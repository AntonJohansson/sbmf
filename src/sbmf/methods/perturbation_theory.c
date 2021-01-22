/*
 * Holds all information needed to do the PT,
 * easy to pass around. Basicly params to
 * pt_rayleigh_schroedinger
 */
struct pt_settings {
	struct nlse_result* res;
	f64* g0;
	i64* particle_count;

	struct eigen_result_real* states;
	const u32 N; /* states to include */
	const u32 L; /* coeff count */
	nlse_operator_func* pert;
};

static inline f64 G0(struct pt_settings* pt, u32 A, u32 B) {
	return pt->g0[A*pt->res->component_count + B];
}

static inline f64 E(struct pt_settings* pt, u32 A, u32 i) {
	return pt->states[A].eigenvalues[i];
}

static inline f64* PHI(struct pt_settings* pt, u32 A, u32 i) {
	return &pt->states[A].eigenvectors[i * pt->L];
}

struct V_params {
	u32 coeff_count;
	f64* i;
	f64* j;
	f64* k;
	f64* l;
};

void V_integrand(f64* out, f64* in, u32 len, void* data) {
	struct V_params* p = data;

	f64 sample_i[len]; ho_sample(p->coeff_count, p->i, len, sample_i, in);
	f64 sample_j[len]; ho_sample(p->coeff_count, p->j, len, sample_j, in);
	f64 sample_k[len]; ho_sample(p->coeff_count, p->k, len, sample_k, in);
	f64 sample_l[len]; ho_sample(p->coeff_count, p->l, len, sample_l, in);

	for (u32 i = 0; i < len; ++i) {
		out[i] = sample_i[i]*sample_j[i]*sample_k[i]*sample_l[i];
	}
}

static inline f64 V(struct pt_settings* pt, u32 A, u32 B, u32 i, u32 j, u32 k, u32 l) {
	/* Check if we have computed this integral already */
	struct V_params p = {
		.coeff_count = pt->L,
		.i = PHI(pt, A, i),
		.j = PHI(pt, B, j),
		.k = PHI(pt, A, k),
		.l = PHI(pt, B, l)
	};

    struct quadgk_settings settings = {
        .abs_error_tol = 1e-15,
        .rel_error_tol = 1e-8,
		.gk = gk15,
        .max_evals = 1e5,
		.userdata = &p,
    };

	struct quadgk_result res;
	quadgk_infinite_interval(V_integrand, &settings, &res);
	assert(res.converged);

	return res.integral;
}

/*
 * Helper functions since these will be calculated a lot
 */

static inline f64 rs_2nd_order_me(struct pt_settings* pt, u32 A, u32 B, u32 i, u32 j) {
	f64 me = 0.0;
	if (A == B) {
		f64 factor = (i == j) ? 1.0/sqrt(2.0) : 1.0;
		me = factor * G0(pt,A,A) * sqrt(pt->particle_count[A] * (pt->particle_count[A] - 1))
			* V(pt, A,A, i,j,0,0);
	} else {
		me = G0(pt,A,B) * sqrt(pt->particle_count[A] * pt->particle_count[B])
			* V(pt, A,B, i,j,0,0);
	}
	return me;
}

static inline f64 rs_2nd_order_ediff(struct pt_settings* pt, u32 A, u32 B, u32 i, u32 j) {
	return E(pt,A,0) + E(pt,B,0) - E(pt,A,i) - E(pt,B,j);
}

/*
 * Main function for Rayleigh-Schrodinger perturbation theory
 */

struct pt_result rayleigh_schroedinger_pt_rf(struct nlse_settings settings, struct nlse_result res, u32 component, f64* g0, i64* particle_count) {
	/* order of hamiltonians, that is include all states */
	const u32 states_to_include = res.coeff_count;
	const i64* N = particle_count;
	sbmf_log_info("running 1comp RSPT:\n    components: %u\n    states: %u\n", res.component_count, states_to_include);

	struct eigen_result_real states[res.component_count];
	for (u32 i = 0; i < res.component_count; ++i) {
		states[i] = find_eigenpairs_full_real(res.hamiltonian[i]);
		//states[i] = find_eigenpairs_sparse_real(res.hamiltonian[i], states_to_include, EV_SMALLEST);
		for (u32 j = 0; j < states_to_include; ++j) {
			f64_normalize(&states[i].eigenvectors[j*res.coeff_count], &states[i].eigenvectors[j*res.coeff_count], res.coeff_count);
		}
	}

	struct pt_settings pt = {
		.res = &res,
		.g0 = g0,
		.particle_count = particle_count,
		.states = states,
		.N = states_to_include,
		.L = res.coeff_count,
		.pert = settings.spatial_pot_perturbation,
	};

	/* Zeroth order PT */
	sbmf_log_info("Starting zeroth order PT");
	f64 E0 = 0.0;
	{
		E0 += N[component] * E(&pt, component, 0);
	}
	sbmf_log_info("\tE0: %e", E0);

	/* First order PT */
	sbmf_log_info("Starting first order PT");
	f64 E1 = 0.0;
	{
		/* Handles interaction within component */
		E1 += -0.5 * G0(&pt,component,component) * N[component] * (N[component]-1) * V(&pt, component,component, 0,0,0,0);
	}
	sbmf_log_info("\tE1: %e", E1);

	const u32 pt2_cache_size = ((states_to_include-1)*(states_to_include))/2;
	f64 pt2_cache[pt2_cache_size];

	/* Assumes i in [0,states_to_include), j in [0,states_to_include) */
#define PT2_CACHE_INDEX(i, j) \
	((i)*states_to_include - (((i)*(i+1))/2) + j)


	/* Second order PT */
	sbmf_log_info("Starting second order PT");
	f64 E2 = 0.0;
	{
		/*
		 * Double substitutions (both excitations within same component),
		 * loop over unique pairs (j,k).
		 */
#pragma omp parallel for reduction(+: E2)
		for (u32 i = 1; i < states_to_include; ++i) {
			for (u32 j = i; j < states_to_include; ++j) {
				f64 me = rs_2nd_order_me(&pt, component,component, i,j);
				f64 Ediff = rs_2nd_order_ediff(&pt, component,component, i,j);

				pt2_cache[PT2_CACHE_INDEX(i-1, j-i)] = me/Ediff;

				E2 += me*me/(Ediff);
			}
		}
	}
	sbmf_log_info("\tE2: %e", E2);

	/* Third order PT */
	sbmf_log_info("Starting third order PT");
	f64 E3 = 0.0;
	{
		f64 E_00_00 = 0;
		{
#pragma omp parallel for reduction(+: E_00_00)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 n = m; n < states_to_include; ++n) {
					f64 tmn = pt2_cache[PT2_CACHE_INDEX(m-1,n-m)];
					E_00_00 += tmn*tmn;
				}
			}

			f64 v_00_00 = G0(&pt,component,component)*V(&pt, component,component, 0,0,0,0);
			E_00_00 *= v_00_00;
		}
		sbmf_log_info("\t\t00,00: %.10e", E_00_00);

		f64 E_m0_n0 = 0;
		{
			const f64 c_root_2_minus_1 = sqrt(2.0) - 1.0;
			const f64 c_3_minus_2_root_2 = 3.0 - 2.0*sqrt(2.0);

#pragma omp parallel for reduction(+: E_m0_n0)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 n = 1; n < states_to_include; ++n) {
					const f64 v_m0_n0 = G0(&pt,component,component)*V(&pt, component,component, m,0,n,0);

					f64 sum = 0.0;
					for (u32 p = 1; p < states_to_include; ++p) {

						f64 tmp = 0;
						if (p >= m)
							tmp = pt2_cache[PT2_CACHE_INDEX(m-1, p-m)];
						else
							tmp = pt2_cache[PT2_CACHE_INDEX(p-1, m-p)];

						f64 tnp = 0;
						if (p >= n)
							tnp = pt2_cache[PT2_CACHE_INDEX(n-1, p-n)];
						else
							tnp = pt2_cache[PT2_CACHE_INDEX(p-1, n-p)];

						const f64 delta_mp = (m == p) ? 1.0 : 0.0;
						const f64 delta_np = (n == p) ? 1.0 : 0.0;

						const f64 coeff = 1 + c_root_2_minus_1*(delta_mp + delta_np) + c_3_minus_2_root_2*(delta_mp*delta_np);
						sum += coeff * tmp * tnp;
					}

					E_m0_n0 += v_m0_n0 * sum;
				}
			}

			E_m0_n0 *= (N[component] - 3);
		}
		sbmf_log_info("\t\tm0,n0: %.10e", E_m0_n0);

		f64 E_mn_pq = 0;
		{
			const f64 c_root_2_minus_2 = sqrt(2.0) - 2.0;
			const f64 c_3_minus_2_root_2 = 3.0 - 2.0*sqrt(2.0);

#pragma omp parallel for reduction(+: E_mn_pq)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 n = m; n < states_to_include; ++n) {
					//if (((m ^ n) & 1) != 0)
					//	continue;

					f64 tmn = pt2_cache[PT2_CACHE_INDEX(m-1,n-m)];
					const f64 delta_mn = (m == n) ? 1.0 : 0.0;

					for (u32 p = 1; p < states_to_include; ++p) {
						for (u32 q = p; q < states_to_include; ++q) {
							//if (((m ^ n) & 1) != 0)
							//	continue;

							f64 v_mn_pq = G0(&pt,component,component)*V(&pt, component,component, m,n,p,q);
							f64 tpq = pt2_cache[PT2_CACHE_INDEX(p-1,q-p)];

							const f64 delta_pq = (p == q) ? 1.0 : 0.0;
							const f64 coeff = 2.0 + c_root_2_minus_2*(delta_mn + delta_pq) + c_3_minus_2_root_2*(delta_mn*delta_pq);
							E_mn_pq += coeff*tmn*tpq*v_mn_pq;
						}
					}

				}
			}
		}
		sbmf_log_info("\t\tmn,pq: %.10e", E_mn_pq);

		E3 = E_00_00 + E_m0_n0 + E_mn_pq;
	}
	sbmf_log_info("\tE3: %e", E3);

	return (struct pt_result) {
		.E0 = E0,
		.E1 = E1,
		.E2 = E2,
		.E3 = E3,
	};
}

struct pt_result rayleigh_schroedinger_pt_rf_2comp(struct nlse_settings settings, struct nlse_result res, f64* g0, i64* particle_count) {
	/* order of hamiltonians, that is include all states */
	const u32 states_to_include = res.coeff_count;
	const i64* N = particle_count;
	sbmf_log_info("running 2comp RSPT:\n    components: %u\n    states: %u\n", res.component_count, states_to_include);

	struct eigen_result_real states[res.component_count];
	for (u32 i = 0; i < res.component_count; ++i) {
		states[i] = find_eigenpairs_full_real(res.hamiltonian[i]);
		for (u32 j = 0; j < states_to_include; ++j) {
			f64_normalize(&states[i].eigenvectors[j*res.coeff_count], &states[i].eigenvectors[j*res.coeff_count], res.coeff_count);
		}
	}

	struct pt_settings pt = {
		.res = &res,
		.g0 = g0,
		.particle_count = particle_count,
		.states = states,
		.N = states_to_include,
		.L = res.coeff_count,
		.pert = settings.spatial_pot_perturbation,
	};

	struct pt_result res_comp_A = rayleigh_schroedinger_pt_rf(settings, res, 0, g0, particle_count);
	struct pt_result res_comp_B = rayleigh_schroedinger_pt_rf(settings, res, 1, g0, particle_count);

	/* Zeroth order PT */
	sbmf_log_info("Starting zeroth order PT");
	f64 E0 = 0.0;
	{
		E0 = res_comp_A.E0 + res_comp_B.E0;
	}
	sbmf_log_info("\tE0: %e", E0);

	/* First order PT */
	sbmf_log_info("Starting first order PT");
	f64 E1 = 0.0;
	{
		E1 = res_comp_A.E1 + res_comp_B.E1;
		/* Handles interaction between components */
		for (u32 A = 0; A < res.component_count; ++A) {
			for (u32 B = A+1; B < res.component_count; ++B) {
				E1 += -G0(&pt,A,B) * N[A] * N[B] * V(&pt, A,B, 0,0,0,0);
			}
		}
	}
	sbmf_log_info("\tE1: %e", E1);

	const u32 pt2_cache_size = (states_to_include-1)*(states_to_include-1);
	f64 pt2_cache_AB[pt2_cache_size];

	/* Assumes i in [0,states_to_include), j in [0,states_to_include) */
#define PT2_CACHE_AB_INDEX(i, j) \
	(i)*(states_to_include-1) + (j)

	/* Second order PT */
	sbmf_log_info("Starting second order PT");
	f64 E2 = 0.0;
	{
		E2 = res_comp_A.E2 + res_comp_B.E2;
		/*
		 * Double substitutions (separate components).
		 * A,B refers to components.
		 */
		for (u32 A = 0; A < res.component_count; ++A) {
			for (u32 B = A+1; B < res.component_count; ++B) {
				/*
				 * loop over all combinations of states in components
				 * A,B.
				 */
#pragma omp parallel for reduction(+: E2)
				for (u32 i = 1; i < states_to_include; ++i) {
					for (u32 j = 1; j < states_to_include; ++j) {
						f64 me = rs_2nd_order_me(&pt, A,B, i,j);
						f64 Ediff = rs_2nd_order_ediff(&pt, A,B, i,j);

						pt2_cache_AB[PT2_CACHE_AB_INDEX(i-1, j-1)] = me/Ediff;

						E2 += me*me/(Ediff);
					}
				}
			}
		}
	}
	sbmf_log_info("\tE2: %e", E2);

	sbmf_log_info("Starting third order PT");
	f64 E3 = 0.0;
	{
		E3 = res_comp_A.E3 + res_comp_B.E3;

		f64 E_00_00 = 0.0;
		{

			f64 sum = 0;
#pragma omp parallel for reduction(+: sum)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 n = 1; n < states_to_include; ++n) {
					const f64 tmn = pt2_cache_AB[PT2_CACHE_AB_INDEX(m-1, n-1)];
					sum += tmn*tmn;
				}
			}

			for (u32 A = 0; A < res.component_count; ++A) {
				for (u32 B = A+1; B < res.component_count; ++B) {
					f64 v_00_00 = V(&pt, A,B, 0,0,0,0);
					E_00_00 += G0(&pt, A,B) * v_00_00 * sum;
				}
			}

		}
		sbmf_log_info("\t\t00,00: %.10e", E_00_00);

		f64 E_m0_n0 = 0.0;
		{
			const u32 A = 0;
			const u32 B = 1;

			{
#pragma omp parallel for reduction(+: E_m0_n0)
				for (u32 m = 1; m < states_to_include; ++m) {
					for (u32 n = 1; n < states_to_include; ++n) {

						f64 sum = 0;
						for (u32 p = 1; p < states_to_include; ++p) {
							const f64 tmp = pt2_cache_AB[PT2_CACHE_AB_INDEX(m-1, p-1)];
							const f64 tnp = pt2_cache_AB[PT2_CACHE_AB_INDEX(n-1, p-1)];
							sum += tmp*tnp;
						}

						const f64 v_AA_m0_n0 = G0(&pt, A,A) * (N[A]-1) * V(&pt, A,A, m,0,n,0);
						const f64 v_AB_m0_n0 = - G0(&pt, A,B) * V(&pt, A,B, m,0,n,0);
						E_m0_n0 += (v_AA_m0_n0 + v_AB_m0_n0) * sum;
					}
				}
			}

			{
#pragma omp parallel for reduction(+: E_m0_n0)
				for (u32 m = 1; m < states_to_include; ++m) {
					for (u32 n = 1; n < states_to_include; ++n) {

						f64 sum = 0;
						for (u32 p = 1; p < states_to_include; ++p) {
							const f64 tpm = pt2_cache_AB[PT2_CACHE_AB_INDEX(p-1, m-1)];
							const f64 tpn = pt2_cache_AB[PT2_CACHE_AB_INDEX(p-1, n-1)];
							sum += tpm*tpn;
						}

						const f64 v_BB_m0_n0 = G0(&pt, B,B) * (N[B]-1) * V(&pt, B,B, m,0,n,0);
						const f64 v_AB_0m_0n = - G0(&pt, A,B) * V(&pt, A,B, 0,m,0,n);
						E_m0_n0 += (v_BB_m0_n0 + v_AB_0m_0n) * sum;
					}
				}
			}
		}
		sbmf_log_info("\t\tm0,n0: %.10e", E_m0_n0);

		f64 E_mn_pq = 0;
		{
			const u32 A = 0;
			const u32 B = 1;

#pragma omp parallel for reduction(+: E_mn_pq)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 n = 1; n < states_to_include; ++n) {
					for (u32 p = 1; p < states_to_include; ++p) {
						for (u32 q = 1; q < states_to_include; ++q) {
							const f64 tmn = pt2_cache_AB[PT2_CACHE_AB_INDEX(m-1, n-1)];
							const f64 tpq = pt2_cache_AB[PT2_CACHE_AB_INDEX(p-1, q-1)];
							const f64 v_mn_pq = G0(&pt, A,B) * V(&pt, A,B, m,n,p,q);
							E_mn_pq += v_mn_pq*tmn*tpq;
						}
					}
				}
			}
		}
		sbmf_log_info("\t\tmn,pq: %.10e", E_mn_pq);

		E3 += E_00_00 + E_m0_n0 + E_mn_pq;
	}
	sbmf_log_info("\tE3: %e", E3);

	return (struct pt_result) {
		.E0 = E0,
		.E1 = E1,
		.E2 = E2,
		.E3 = E3,
	};
}















struct Vp_params {
	u32 coeff_count;
	f64* i;
	f64* j;
	nlse_operator_func* pert;
};

void Vp_integrand(f64* out, f64* in, u32 len, void* data) {
	struct Vp_params* p = data;

	f64 sample_i[len]; ho_sample(p->coeff_count, p->i, len, sample_i, in);
	f64 sample_j[len]; ho_sample(p->coeff_count, p->j, len, sample_j, in);

	f64 sample_pert[len];
	p->pert(len, sample_pert, in, 0, NULL, NULL);

	for (u32 i = 0; i < len; ++i) {
		out[i] = sample_i[i]*sample_pert[i]*sample_j[i];
	}
}

static inline f64 Vp(struct pt_settings* pt, u32 A, u32 i, u32 j) {
	/* Check if we have computed this integral already */
	struct Vp_params p = {
		.coeff_count = pt->L,
		.i = PHI(pt, A, i),
		.j = PHI(pt, A, j),
		.pert = pt->pert,
	};

    struct quadgk_settings settings = {
        .abs_error_tol = 1e-15,
        .rel_error_tol = 1e-8,
		.gk = gk15,
        .max_evals = 1e5,
		.userdata = &p,
    };

	struct quadgk_result res;
	quadgk_infinite_interval(Vp_integrand, &settings, &res);
	assert(res.converged);

	return res.integral;
}



static inline f64 en_nhn(struct pt_settings* pt, u32 A, u32 i, u32 j) {
	f64 sum = 0;
	for (u32 k = 0; k < pt->L; ++k) {
		sum += pt->states[A].eigenvectors[i*pt->L + k] * pt->states[A].eigenvectors[j*pt->L + k] * ho_eigenval(k);

	}
	/* <i|Vp|j> */
	if (pt->pert)
		sum += Vp(pt, A, i, j);

	return sum;
}



static inline f64 en_nHn(struct pt_settings* pt, u32 A, u32 i, u32 j) {
	if (i == j) {
		if (i == 0) {
			/* <00|H|00> */
			return pt->particle_count[A] * en_nhn(pt, A, 0, 0)
				+ 0.5 * G0(pt,A,A) * pt->particle_count[A]*(pt->particle_count[A]-1) * V(pt, A,A, 0,0,0,0);
		}

		/* <mm|H|mm> */
		return 2*en_nhn(pt, A, i,i) + (pt->particle_count[A]-2) * en_nhn(pt, A, 0,0)
			+ 0.5 * G0(pt,A,A) * (
						2*V(pt, A,A, i,i,i,i)
						+8*(pt->particle_count[A]-2)*V(pt, A,A, i,0,i,0)
						+(pt->particle_count[A]-2)*(pt->particle_count[A]-3)*V(pt, A,A, 0,0,0,0)
					);
	} else {
		/* <mn|H|mn> */
		return en_nhn(pt, A, i,i) + en_nhn(pt, A, j,j) + (pt->particle_count[A]-2) * en_nhn(pt, A, 0,0)
			+ 0.5 * G0(pt,A,A) * (
						4*V(pt, A,A, i,j,i,j)
						+4*(pt->particle_count[A]-2)*V(pt, A,A, i,0,i,0)
						+4*(pt->particle_count[A]-2)*V(pt, A,A, j,0,j,0)
						+(pt->particle_count[A]-2)*(pt->particle_count[A]-3)*V(pt, A,A, 0,0,0,0)
					);
	}
}










struct pt_result en_pt_rf(struct nlse_settings settings, struct nlse_result res, u32 component, f64* g0, i64* particle_count) {
	/* order of hamiltonians, that is include all states */
	const u32 states_to_include = res.coeff_count;
	const i64* N = particle_count;
	sbmf_log_info("running 1comp ENPT:\n    components: %u\n    states: %u\n", res.component_count, states_to_include);

	struct eigen_result_real states[res.component_count];
	for (u32 i = 0; i < res.component_count; ++i) {
		states[i] = find_eigenpairs_full_real(res.hamiltonian[i]);
		//states[i] = find_eigenpairs_sparse_real(res.hamiltonian[i], states_to_include, EV_SMALLEST);
		for (u32 j = 0; j < states_to_include; ++j) {
			f64_normalize(&states[i].eigenvectors[j*res.coeff_count], &states[i].eigenvectors[j*res.coeff_count], res.coeff_count);
		}
	}


	struct pt_settings pt = {
		.res = &res,
		.g0 = g0,
		.particle_count = particle_count,
		.states = states,
		.N = states_to_include,
		.L = res.coeff_count,
		.pert = settings.spatial_pot_perturbation,
	};

	/* Zeroth order PT */
	sbmf_log_info("Starting zeroth order PT");
	f64 E0 = 0.0;
	{
		E0 += en_nHn(&pt, component, 0,0);
	}
	sbmf_log_info("\tE0: %e", E0);

	/* E1 always zero in EN PT */
	f64 E1 = 0.0;

	const u32 pt2_cache_size = ((states_to_include-1)*(states_to_include))/2;
	f64 pt2_cache[pt2_cache_size];

	/* Assumes i in [0,states_to_include), j in [0,states_to_include) */
#define PT2_CACHE_INDEX(i, j) \
	((i)*states_to_include - (((i)*(i+1))/2) + j)


	/* Second order PT */
	sbmf_log_info("Starting second order PT");
	f64 E2 = 0.0;
	{
		/*
		 * Double substitutions (both excitations within same component),
		 * loop over unique pairs (j,k).
		 */
		const f64 E_0H0 = en_nHn(&pt, component, 0,0);
#pragma omp parallel for reduction(+: E2)
		for (u32 i = 1; i < states_to_include; ++i) {
			for (u32 j = i; j < states_to_include; ++j) {
				f64 me = rs_2nd_order_me(&pt, component,component, i,j);
				f64 Ediff = E_0H0 - en_nHn(&pt, component, i,j);
				//sbmf_log_info("-------: %lf --- %lf", me*me, Ediff);

				pt2_cache[PT2_CACHE_INDEX(i-1, j-i)] = me/Ediff;

				E2 += me*me/(Ediff);
			}
		}
	}
	sbmf_log_info("\tE2: %e", E2);

	/* Third order PT */
	sbmf_log_info("Starting third order PT");
	f64 E3 = 0.0;
	{
		/* Comes from terms of the form <k|V|k> which are always zero in EN PT */
		f64 E_00_00 = 0;

		f64 E_m0_n0 = 0;
		{
			const f64 c_root_2_minus_1 = sqrt(2.0) - 1.0;
			const f64 c_3_minus_2_root_2 = 3.0 - 2.0*sqrt(2.0);

#pragma omp parallel for reduction(+: E_m0_n0)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 n = 1; n < states_to_include; ++n) {
					/* avoid <mn|V|mn> = 0 */
					if (m == n)
						continue;

					const f64 v_m0_n0 = G0(&pt,component,component)*V(&pt, component,component, m,0,n,0);

					f64 sum = 0.0;
					for (u32 p = 1; p < states_to_include; ++p) {

						f64 tmp = 0;
						if (p >= m)
							tmp = pt2_cache[PT2_CACHE_INDEX(m-1, p-m)];
						else
							tmp = pt2_cache[PT2_CACHE_INDEX(p-1, m-p)];

						f64 tnp = 0;
						if (p >= n)
							tnp = pt2_cache[PT2_CACHE_INDEX(n-1, p-n)];
						else
							tnp = pt2_cache[PT2_CACHE_INDEX(p-1, n-p)];

						const f64 delta_mp = (m == p) ? 1.0 : 0.0;
						const f64 delta_np = (n == p) ? 1.0 : 0.0;

						const f64 coeff = 1 + c_root_2_minus_1*(delta_mp + delta_np) + c_3_minus_2_root_2*(delta_mp*delta_np);
						sum += coeff * tmp * tnp;
					}

					E_m0_n0 += v_m0_n0 * sum;

				}
			}

			E_m0_n0 *= (N[component] - 3);
		}
		sbmf_log_info("\t\tm0,n0: %.10e", E_m0_n0);

		f64 E_mn_pq = 0;
		{
			const f64 c_root_2_minus_2 = sqrt(2.0) - 2.0;
			const f64 c_3_minus_2_root_2 = 3.0 - 2.0*sqrt(2.0);

#pragma omp parallel for reduction(+: E_mn_pq)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 n = m; n < states_to_include; ++n) {
					if (((m ^ n) & 1) != 0)
						continue;

					f64 tmn = pt2_cache[PT2_CACHE_INDEX(m-1,n-m)];
					const f64 delta_mn = (m == n) ? 1.0 : 0.0;

					for (u32 p = 1; p < states_to_include; ++p) {
						for (u32 q = p; q < states_to_include; ++q) {
							if (((m ^ n) & 1) != 0)
								continue;
							if (m*states_to_include + n == p*states_to_include + q)
								continue;

							f64 v_mn_pq = G0(&pt,component,component)*V(&pt, component,component, m,n,p,q);
							f64 tpq = pt2_cache[PT2_CACHE_INDEX(p-1,q-p)];

							const f64 delta_pq = (p == q) ? 1.0 : 0.0;
							const f64 coeff = 2.0 + c_root_2_minus_2*(delta_mn + delta_pq) + c_3_minus_2_root_2*(delta_mn*delta_pq);
							E_mn_pq += coeff*tmn*tpq*v_mn_pq;
						}
					}

				}
			}
		}
		sbmf_log_info("\t\tmn,pq: %.10e", E_mn_pq);

		E3 = E_00_00 + E_m0_n0 + E_mn_pq;
	}
	sbmf_log_info("\tE3: %e", E3);

	return (struct pt_result) {
		.E0 = E0,
		.E1 = E1,
		.E2 = E2,
		.E3 = E3,
	};
}
