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

static f64 V(const u32 coeff_count, f64 i[static coeff_count],
							 f64 j[static coeff_count],
							 f64 k[static coeff_count],
							 f64 l[static coeff_count]) {
	struct V_params p = {
		.coeff_count = coeff_count,
		.i = i,
		.j = j,
		.k = k,
		.l = l
	};

    struct integration_settings settings = {
        .abs_error_tol = 1e-12,
        .rel_error_tol = 1e-10,
		.gk = gk15,
        .max_evals = 1e5,
		.userdata = &p,
    };

	struct integration_result res = quadgk(V_integrand, -INFINITY, INFINITY, settings);
	assert(res.converged);

	return res.integral;
}

struct pt_result rayleigh_schroedinger_pt(struct nlse_result res, f64* g0, u32* particle_count) {
	/* order of hamiltonians, that is include all states */
	const u32 states_to_include = res.coeff_count;
	sbmf_log_info("running rayleigh schödinger PT:\n    components: %u\n    states: %u\n", res.component_count, states_to_include);

	struct eigen_result_real states[res.component_count];
	for (u32 i = 0; i < res.component_count; ++i) {
		//states[i] = find_eigenpairs_sparse_real(res.hamiltonian[i], states_to_include, EV_SMALLEST_MAG);
		states[i] = find_eigenpairs_full_real(res.hamiltonian[i]);
		for (u32 j = 0; j < states_to_include; ++j) {
			f64_normalize(&states[i].eigenvectors[j*res.coeff_count], &states[i].eigenvectors[j*res.coeff_count], res.coeff_count);
		}
	}

	const u32 L = res.coeff_count;

	/*
	 * Macros are not local, but they are only used here...
	 * PHI returns the coeffs. of state 'state' in component
	 * 'component'. ENERGY returns the energy of that same
	 * state. G0 return the interaction strength of component
	 * A with respect to component B, the g0 array is assumed
	 * to be symmetric.
	 */
#define PHI(component, state) 	 states[component].eigenvectors[state * res.coeff_count]
#define ENERGY(component, state) states[component].eigenvalues[state]
#define G0(componentA, componentB) 	 g0[componentA*res.component_count + componentB]

	/* zeroth order PT, computes <0|H0|0> */
	f64 E0 = 0.0;
	{
		for (u32 i = 0; i < res.component_count; ++i) {
			E0 += particle_count[i] * ENERGY(i,0);
		}
	}

	/* first ordet PT, computes <0|V|0> */
	f64 E1 = 0.0;
	{
		/* Handles interaction within component */
		for (u32 i = 0; i < res.component_count; ++i) {
			E1 += -0.5 * G0(i,i) * particle_count[i] * (particle_count[i]-1)
				* V(L, &PHI(i,0), &PHI(i,0), &PHI(i,0), &PHI(i,0));
		}

		/* Handles interaction between components */
		for (u32 i = 0; i < res.component_count; ++i) {
			for (u32 j = i+1; j < res.component_count; ++j) {
				E1 += -G0(i,j) * particle_count[i] * particle_count[j]
					* V(L, &PHI(i,0), &PHI(j,0), &PHI(i,0), &PHI(j,0));

			}
		}
	}

	/* second order pt, computes
	 *     sum |<k|V|0>|^2/(E0-Ek)
	 * where |k> are double substitions.
	 */
	f64 E2 = 0.0;
	f64 E3_last_term = 0.0;
	{
		for (u32 A = 0; A < res.component_count; ++A) {
			/* Double substitutions (both excitations within same component),
			 * loop over unique pairs (j,k).
			 */
			for (u32 i = 1; i < states_to_include; ++i) {
				for (u32 j = i; j < states_to_include; ++j) {
					f64 factor = (i == j) ? 1.0/sqrt(2.0) : 1.0;
					f64 me = factor * G0(A,A) * sqrt(particle_count[A] * (particle_count[A] - 1))
						* V(L, &PHI(A,i), &PHI(A,j), &PHI(A,0), &PHI(A,0));
					f64 Ediff = 2*ENERGY(A,0) - (ENERGY(A,i) + ENERGY(A,j));
					E2 += me*me/(Ediff);
					E3_last_term += me*me/(Ediff*Ediff);
				}
			}
		}

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
				for (u32 i = 1; i < states_to_include; ++i) {
					for (u32 j = 1; j < states_to_include; ++j) {
						f64 me = G0(A,B) * sqrt(particle_count[A] * particle_count[B])
								* V(L, &PHI(A,i), &PHI(B,j), &PHI(A,0), &PHI(B,0));
						f64 Ediff = ENERGY(A,0) + ENERGY(B,0) - (ENERGY(A,i) + ENERGY(B,j));
						E2 += me*me/(Ediff);
						E3_last_term += me*me/(Ediff*Ediff);
					}
				}
			}
		}
	}

	f64 E3 = 0.0;
#if 0
	{
		/*
		 * Loop over all <AmAn|V|ApAq>,
		 * requires all unique combinations (m,n,p,q)
		 */

		/* mm,mm */
		sbmf_log_info("mm,mm");
		for (u32 m = 1; m < states_to_include; ++m) {
			//sbmf_log_info("(%u,%u),(%u,%u)", m,m,m,m);
			f64 v_AA_mmmm = V(L, &PHI(0,m), &PHI(0,m), &PHI(0,m), &PHI(0,m));
			f64 v_AA_m0m0 = V(L, &PHI(0,m), &PHI(0,0), &PHI(0,m), &PHI(0,0));
			f64 v_AA_0000 = V(L, &PHI(0,0), &PHI(0,0), &PHI(0,0), &PHI(0,0));
			f64 v_AA_mm00 = V(L, &PHI(0,m), &PHI(0,m), &PHI(0,0), &PHI(0,0));
			//f64 v_BA_0000 = V(L, &PHI(1,0), &PHI(0,0), &PHI(1,0), &PHI(0,0));

			/* Two double substitutions */
			f64 me0 = 0.5*G0(0,0)*(2*v_AA_mmmm
						+ 8*(particle_count[0]-2)*v_AA_m0m0
						+ (particle_count[0]-2)*(particle_count[0]-3)*v_AA_0000)
					- G0(0,0)*(particle_count[0]-1)*(2*v_AA_m0m0 + (particle_count[0]-2)*v_AA_0000);
					//- G0(0,1)*particle_count[0]*particle_count[1]*v_BA_0000;

			/* One double substitution */
			f64 me1;
			{
				f64 factor = (m == m) ? 1.0/sqrt(2.0) : 1.0;
				me1 = factor * G0(0,0) * sqrt(particle_count[0] * (particle_count[0] - 1))
					* v_AA_mm00;
			}

			f64 Ediff = 2*ENERGY(0,0) - 2*ENERGY(0,m);
			E3 += (me1*me0*me1)/(Ediff*Ediff);
		}

		/* mn,mn */
		sbmf_log_info("mn,mn");
		for (u32 m = 1; m < states_to_include; ++m) {
			for (u32 n = m+1; n < states_to_include; ++n) {
				//sbmf_log_info("(%u,%u),(%u,%u)", m,n,m,n);
				f64 v_AA_mnmn = V(L, &PHI(0,m), &PHI(0,n), &PHI(0,m), &PHI(0,n));
				f64 v_AA_m0m0 = V(L, &PHI(0,m), &PHI(0,0), &PHI(0,m), &PHI(0,0));
				f64 v_AA_n0n0 = V(L, &PHI(0,n), &PHI(0,0), &PHI(0,n), &PHI(0,0));
				f64 v_AA_0000 = V(L, &PHI(0,0), &PHI(0,0), &PHI(0,0), &PHI(0,0));
				f64 v_AA_mn00 = V(L, &PHI(0,m), &PHI(0,n), &PHI(0,0), &PHI(0,0));
				//f64 v_BA_0000 = V(L, &PHI(1,0), &PHI(0,0), &PHI(1,0), &PHI(0,0));

				/* Two double substitutions */
				f64 me0 = 0.5*G0(0,0)*(
						4*v_AA_mnmn
						+ 4*(particle_count[0]-2)*v_AA_m0m0
						+ 4*(particle_count[0]-2)*v_AA_n0n0
						+ (particle_count[0]-2)*(particle_count[0]-3)*v_AA_0000)
					- G0(0,0)*(particle_count[0]-1)*(v_AA_m0m0 + v_AA_n0n0 + (particle_count[0]-2)*v_AA_0000);
				//- G0(0,1)*particle_count[0]*particle_count[1]*v_BA_0000;

				/* One double substitution */
				f64 me1;
				{
					f64 factor = (m == n) ? 1.0/sqrt(2.0) : 1.0;
					me1 = factor * G0(0,0) * sqrt(particle_count[0] * (particle_count[0] - 1))
						* v_AA_mn00;
				}

				f64 Ediff = 2*ENERGY(0,0) - ENERGY(0,m) - ENERGY(0,n);
				E3 += (me1*me0*me1)/(Ediff*Ediff);
			}
		}

		/* mm,nn */
		sbmf_log_info("mm,nn");
		for (u32 m = 1; m < states_to_include; ++m) {
			for (u32 n = m+1; n < states_to_include; ++n) {
				//sbmf_log_info("(%u,%u),(%u,%u)", m,m,n,n);
				f64 v_AA_mmnn = V(L, &PHI(0,m), &PHI(0,m), &PHI(0,n), &PHI(0,n));
				f64 v_AA_mm00 = V(L, &PHI(0,m), &PHI(0,m), &PHI(0,0), &PHI(0,0));
				f64 v_AA_nn00 = V(L, &PHI(0,n), &PHI(0,n), &PHI(0,0), &PHI(0,0));

				/* Two double substitutions */
				f64 me0 = G0(0,0) * v_AA_mmnn;
				/* One double substitution */
				f64 me1;
				{
					f64 factor = (m == m) ? 1.0/sqrt(2.0) : 1.0;
					me1 = factor * G0(0,0) * sqrt(particle_count[0] * (particle_count[0] - 1))
						* v_AA_mm00;
				}
				f64 me2;
				{
					f64 factor = (n == n) ? 1.0/sqrt(2.0) : 1.0;
					me2 = factor * G0(0,0) * sqrt(particle_count[0] * (particle_count[0] - 1))
						* v_AA_nn00;
				}

				f64 Ediff0 = 2*ENERGY(0,0) - ENERGY(0,m) - ENERGY(0,m);
				f64 Ediff1 = 2*ENERGY(0,0) - ENERGY(0,n) - ENERGY(0,n);
				E3 += (me1*me0*me2)/(Ediff0*Ediff1);
			}
		}
		/* mm,mp */
		sbmf_log_info("mm,mp");
		for (u32 m = 1; m < states_to_include; ++m) {
			for (u32 p = 1; p < states_to_include; ++p) {
				if (m == p)
					continue;
				//sbmf_log_info("(%u,%u),(%u,%u)", m,m,m,p);
				f64 v_AA_m0p0 = V(L, &PHI(0,m), &PHI(0,0), &PHI(0,p), &PHI(0,0));
				f64 v_AA_mmmp = V(L, &PHI(0,m), &PHI(0,m), &PHI(0,m), &PHI(0,p));
				f64 v_AA_mm00 = V(L, &PHI(0,m), &PHI(0,m), &PHI(0,0), &PHI(0,0));
				f64 v_AA_mp00 = V(L, &PHI(0,m), &PHI(0,p), &PHI(0,0), &PHI(0,0));

				/* Two double substitutions */
				f64 me0 = G0(0,0) * ((particle_count[0]-3)*v_AA_m0p0 + sqrt(2)*v_AA_mmmp);
				/* One double substitution */
				f64 me1;
				{
					f64 factor = (m == m) ? 1.0/sqrt(2.0) : 1.0;
					me1 = factor * G0(0,0) * sqrt(particle_count[0] * (particle_count[0] - 1))
						* v_AA_mm00;
				}
				f64 me2;
				{
					f64 factor = (m == p) ? 1.0/sqrt(2.0) : 1.0;
					me2 = factor * G0(0,0) * sqrt(particle_count[0] * (particle_count[0] - 1))
						* v_AA_mp00;
				}

				f64 Ediff0 = 2*ENERGY(0,0) - ENERGY(0,m) - ENERGY(0,m);
				f64 Ediff1 = 2*ENERGY(0,0) - ENERGY(0,m) - ENERGY(0,p);
				E3 += (me1*me0*me2)/(Ediff0*Ediff1);
			}
		}
		/* mm,pq */
		sbmf_log_info("mm,pq");
		for (u32 m = 1; m < states_to_include; ++m) {
			for (u32 p = 1; p < states_to_include; ++p) {
				if (p == m)
					continue;

				for (u32 q = p+1; q < states_to_include; ++q) {
					if (q == m)
						continue;
					//sbmf_log_info("(%u,%u),(%u,%u)", m,m,p,q);
					f64 v_AA_mmpq = V(L, &PHI(0,m), &PHI(0,m), &PHI(0,p), &PHI(0,q));
					f64 v_AA_mm00 = V(L, &PHI(0,m), &PHI(0,m), &PHI(0,0), &PHI(0,0));
					f64 v_AA_pq00 = V(L, &PHI(0,p), &PHI(0,q), &PHI(0,0), &PHI(0,0));

					/* Two double substitutions */
					f64 me0 = (1.0/sqrt(2.0)) * G0(0,0) * v_AA_mmpq;
					/* One double substitution */
					f64 me1;
					{
						f64 factor = (m == m) ? 1.0/sqrt(2.0) : 1.0;
						me1 = factor * G0(0,0) * sqrt(particle_count[0] * (particle_count[0] - 1))
							* v_AA_mm00;
					}
					f64 me2;
					{
						f64 factor = (p == q) ? 1.0/sqrt(2.0) : 1.0;
						me2 = factor * G0(0,0) * sqrt(particle_count[0] * (particle_count[0] - 1))
							* v_AA_pq00;
					}

					f64 Ediff0 = 2*ENERGY(0,0) - ENERGY(0,m) - ENERGY(0,m);
					f64 Ediff1 = 2*ENERGY(0,0) - ENERGY(0,p) - ENERGY(0,q);
					E3 += (me1*me0*me2)/(Ediff0*Ediff1);
				}
			}
		}

		/* mp,mq */
		sbmf_log_info("mp,mq");
		for (u32 m = 1; m < states_to_include; ++m) {
			for (u32 p = 1; p < states_to_include; ++p) {
				if (p == m)
					continue;
				for (u32 q = p+1; q < states_to_include; ++q) {
					if (q == m)
						continue;
					//sbmf_log_info("(%u,%u),(%u,%u)", m,p,m,q);
					f64 v_AA_p0q0 = V(L, &PHI(0,p), &PHI(0,0), &PHI(0,q), &PHI(0,0));
					f64 v_AA_mpmq = V(L, &PHI(0,m), &PHI(0,p), &PHI(0,m), &PHI(0,q));
					f64 v_AA_mp00 = V(L, &PHI(0,m), &PHI(0,p), &PHI(0,0), &PHI(0,0));
					f64 v_AA_mq00 = V(L, &PHI(0,m), &PHI(0,q), &PHI(0,0), &PHI(0,0));

					/* Two double substitutions */
					f64 me0 = G0(0,0) * ((particle_count[0]-3)*v_AA_p0q0 + 2*v_AA_mpmq);
					/* One double substitution */
					f64 me1;
					{
						f64 factor = (m == p) ? 1.0/sqrt(2.0) : 1.0;
						me1 = factor * G0(0,0) * sqrt(particle_count[0] * (particle_count[0] - 1))
							* v_AA_mp00;
					}
					f64 me2;
					{
						f64 factor = (m == q) ? 1.0/sqrt(2.0) : 1.0;
						me2 = factor * G0(0,0) * sqrt(particle_count[0] * (particle_count[0] - 1))
							* v_AA_mq00;
					}

					f64 Ediff0 = 2*ENERGY(0,0) - ENERGY(0,m) - ENERGY(0,p);
					f64 Ediff1 = 2*ENERGY(0,0) - ENERGY(0,m) - ENERGY(0,q);
					E3 += (me1*me0*me2)/(Ediff0*Ediff1);
				}
			}
		}

		/* mn,pq */
		sbmf_log_info("mn,pq");
		for (u32 m = 1; m < states_to_include; ++m) {
			for (u32 n = m+1; n < states_to_include; ++n) {
				for (u32 p = n+1; p < states_to_include; ++p) {
					for (u32 q = p+1; q < states_to_include; ++q) {
						//sbmf_log_info("(%u,%u),(%u,%u)", m,n,p,q);
						f64 v_AA_mnpq = V(L, &PHI(0,m), &PHI(0,n), &PHI(0,p), &PHI(0,q));
						f64 v_AA_mn00 = V(L, &PHI(0,m), &PHI(0,n), &PHI(0,0), &PHI(0,0));
						f64 v_AA_pq00 = V(L, &PHI(0,p), &PHI(0,q), &PHI(0,0), &PHI(0,0));

						/* Two double substitutions */
						f64 me0 = 2*G0(0,0)*v_AA_mnpq;
						/* One double substitution */
						f64 me1;
						{
							f64 factor = (m == n) ? 1.0/sqrt(2.0) : 1.0;
							me1 = factor * G0(0,0) * sqrt(particle_count[0] * (particle_count[0] - 1))
								* v_AA_mn00;
						}
						f64 me2;
						{
							f64 factor = (p == q) ? 1.0/sqrt(2.0) : 1.0;
							me2 = factor * G0(0,0) * sqrt(particle_count[0] * (particle_count[0] - 1))
								* v_AA_pq00;
						}

						f64 Ediff0 = 2*ENERGY(0,0) - ENERGY(0,m) - ENERGY(0,n);
						f64 Ediff1 = 2*ENERGY(0,0) - ENERGY(0,p) - ENERGY(0,q);
						E3 += (me1*me0*me2)/(Ediff0*Ediff1);
					}
				}
			}
		}

		E3 -= E1 * E3_last_term;
	}
#endif

	return (struct pt_result) {
		.E0 = E0,
		.E1 = E1,
		.E2 = E2,
		.E3 = E3,
	};
}
