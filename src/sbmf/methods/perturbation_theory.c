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

struct pt_result rayleigh_schroedinger_pt(struct nlse_result res, f64* g0, i64* particle_count) {
	/* order of hamiltonians, that is include all states */
	const u32 states_to_include = res.coeff_count;
	sbmf_log_info("running rayleigh schödinger PT:\n    components: %u\n    states: %u\n", res.component_count, states_to_include);

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
	};

	/* zeroth order PT, computes <0|H0|0> */
	f64 E0 = 0.0;
	{
		for (u32 A = 0; A < res.component_count; ++A) {
			E0 += particle_count[A] * E(&pt, A, 0);
		}
	}
	sbmf_log_info("E0: %.10lf", E0);

	/* first ordet PT, computes <0|V|0> */
	f64 E1 = 0.0;
	{
		/* Handles interaction within component */
		for (u32 A = 0; A < res.component_count; ++A) {
			E1 += -0.5 * G0(&pt,A,A) * particle_count[A] * (particle_count[A]-1)
				* V(&pt, A,A, 0,0,0,0);
		}

		/* Handles interaction between components */
		for (u32 A = 0; A < res.component_count; ++A) {
			for (u32 B = A+1; B < res.component_count; ++B) {
				E1 += -G0(&pt,A,B) * particle_count[A] * particle_count[B]
					* V(&pt, A,B, 0,0,0,0);
			}
		}
	}

	/*
	 * second order pt, computes
	 *     sum |<k|V|0>|^2/(E0-Ek)
	 * where |k> are double substitions.
	 */
	f64 E2 = 0.0;
	f64 E3_last_term = 0.0;
	{
		for (u32 A = 0; A < res.component_count; ++A) {
			/*
			 * Double substitutions (both excitations within same component),
			 * loop over unique pairs (j,k).
			 */
			for (u32 i = 1; i < states_to_include; ++i) {
				for (u32 j = i; j < states_to_include; ++j) {
					f64 me = rs_2nd_order_me(&pt, A,A, i,j);
					f64 Ediff = rs_2nd_order_ediff(&pt, A,A, i,j);
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
						f64 me = rs_2nd_order_me(&pt, A,B, i,j);
						f64 Ediff = rs_2nd_order_ediff(&pt, A,B, i,j);
						E2 += me*me/(Ediff);
						E3_last_term += me*me/(Ediff*Ediff);
					}
				}
			}
		}
	}

	sbmf_log_info("\tE0+E1+E2: %.10lf", E0+E1+E2);

	f64 E3 = 0.0;
	f64 E3_before = 0.0;
#if 1
	{
		struct index {u32 m, n, p, q;};
		/*
		 * Loop over all <AmAn|V|ApAq>,
		 * requires all unique combinations (m,n,p,q)
		 */

		/* mm,mm */
		/* states_to_include */
		sbmf_log_info("mm,mm");
		E3_before = E3;
		for (u32 A = 0; A < res.component_count; ++A) {
//#pragma omp parallel for reduction(+: E3)
			for (u32 m = 1; m < states_to_include; ++m) {
				//sbmf_log_info("\t(%u,%u),(%u,%u)", m,m,m,m);
				f64 v_AA_mmmm = V(&pt, A,A, m,m,m,m);
				f64 v_AA_m0m0 = V(&pt, A,A, m,0,m,0);
				f64 v_AA_0000 = V(&pt, A,A, 0,0,0,0);

				/* Two double substitutions */
				f64 me0 = 0.5*G0(&pt,A,A)*(
							2*v_AA_mmmm
							+ 8*(particle_count[A]-2)*v_AA_m0m0
							+ (particle_count[A]-2)*(particle_count[A]-3)*v_AA_0000)
						- G0(&pt,A,A)*(particle_count[A]-1)*(2*v_AA_m0m0 + (particle_count[A]-2)*v_AA_0000);

				/* Handle intracomponent terms */
				for (u32 B = 0; B < res.component_count; ++B) {
					if (B == A)
						continue;

					f64 v_BA_0000 = V(&pt, B,A, 0,0,0,0);
					me0 -= G0(&pt,A,B)*particle_count[A]*particle_count[B]*v_BA_0000;
				}

				/* One double substitution */
				f64 me1 = rs_2nd_order_me(&pt, A,A, m,m);
				f64 Ediff = rs_2nd_order_ediff(&pt, A,A, m,m);

				E3 += (me1*me0*me1)/(Ediff*Ediff);
			}
		}
		sbmf_log_info("\tE3: %lf -- contrib: %lf", E3, E3-E3_before);

		/* mn,mn */
		/* states_to_include * (states_to_include - 1)/2 */
		sbmf_log_info("mn,mn");
		E3_before = E3;
		for (u32 A = 0; A < res.component_count; ++A) {
//#pragma omp parallel for reduction(+: E3)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 n = m+1; n < states_to_include; ++n) {
					//sbmf_log_info("\t(%u,%u),(%u,%u)", m,n,m,n);
					f64 v_AA_mnmn = V(&pt, A,A, m,n,m,n);
					f64 v_AA_m0m0 = V(&pt, A,A, m,0,m,0);
					f64 v_AA_n0n0 = V(&pt, A,A, n,0,n,0);
					f64 v_AA_0000 = V(&pt, A,A, 0,0,0,0);

					/* Two double substitutions */
					f64 me0 = 0.5*G0(&pt,A,A)*(
							4*v_AA_mnmn
							+ 4*(particle_count[A]-2)*v_AA_m0m0
							+ 4*(particle_count[A]-2)*v_AA_n0n0
							+ (particle_count[A]-2)*(particle_count[A]-3)*v_AA_0000)
						- G0(&pt,A,A)*(particle_count[A]-1)*(v_AA_m0m0 + v_AA_n0n0 + (particle_count[A]-2)*v_AA_0000);

					/* Handle intercomponent term */
					for (u32 B = 0; B < res.component_count; ++B) {
						if (B == A)
							continue;

						f64 v_BA_0000 = V(&pt, B,A, 0,0,0,0);
						me0 -= G0(&pt,A,B)*particle_count[A]*particle_count[B]*v_BA_0000;
					}

					/* One double substitution */
					f64 me1 = rs_2nd_order_me(&pt, A,A, m,n);
					f64 Ediff = rs_2nd_order_ediff(&pt, A,A, m,n);
					E3 += (me1*me0*me1)/(Ediff*Ediff);
				}
			}
		}
		sbmf_log_info("\tE3: %lf -- contrib: %lf", E3, E3-E3_before);

		/* mm,nn */
		/* states_to_include * (states_to_include - 1)/2 */
		sbmf_log_info("mm,nn");
		E3_before = E3;
		for (u32 A = 0; A < res.component_count; ++A) {
//#pragma omp parallel for reduction(+: E3)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 n = m+1; n < states_to_include; ++n) {
					//sbmf_log_info("\t(%u,%u),(%u,%u)", m,m,n,n);
					f64 v_AA_mmnn = V(&pt, A,A, m,m,n,n);

					/* Two double substitutions */
					f64 me0 = G0(&pt,A,A) * v_AA_mmnn;

					/* One double substitution */
					f64 me1 = rs_2nd_order_me(&pt, A,A, m,m);
					f64 Ediff0 = rs_2nd_order_ediff(&pt, A,A, m,m);
					f64 me2 = rs_2nd_order_me(&pt, A,A, n,n);
					f64 Ediff1 = rs_2nd_order_ediff(&pt, A,A, n,n);

					E3 += (me1*me0*me2)/(Ediff0*Ediff1);
				}
			}
		}
		sbmf_log_info("\tE3: %lf -- contrib: %lf", E3, E3-E3_before);

		/* mm,mp */
		/* states_to_include * (states_to_include - 1) */
		sbmf_log_info("mm,mp");
		E3_before = E3;
		for (u32 A = 0; A < res.component_count; ++A) {
//#pragma omp parallel for reduction(+: E3) collapse(2)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 p = 1; p < states_to_include; ++p) {
					if (m == p)
						continue;
					//sbmf_log_info("\t(%u,%u),(%u,%u)", m,m,m,p);
					f64 v_AA_m0p0 = V(&pt, A,A, m,0,p,0);
					f64 v_AA_mmmp = V(&pt, A,A, m,m,m,p);

					/* Two double substitutions */
					f64 me0 = G0(&pt,A,A) * ((particle_count[A]-3)*v_AA_m0p0 + sqrt(2)*v_AA_mmmp);

					/* One double substitution */
					f64 me1 	= rs_2nd_order_me(&pt, A,A, m,m);
					f64 Ediff0  = rs_2nd_order_ediff(&pt, A,A, m,m);

					f64 me2 	= rs_2nd_order_me(&pt, A,A, m,p);
					f64 Ediff1  = rs_2nd_order_ediff(&pt, A,A, m,p);

					E3 += (me1*me0*me2)/(Ediff0*Ediff1);
				}
			}
		}
		sbmf_log_info("\tE3: %lf -- contrib: %lf", E3, E3-E3_before);

		/* mm,pq */
		/* (states_to_include-1)*(states_to_include - 2)*(states_to_include-3) */
		sbmf_log_info("mm,pq");
		E3_before = E3;
		{
			const u32 len = (states_to_include-1)*(states_to_include - 2)*(states_to_include-3)/2;
			struct index inds[len];
			u32 iter = 0;
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 p = 1; p < states_to_include; ++p) {
					if (p == m)
						continue;
					for (u32 q = p+1; q < states_to_include; ++q) {
						if (q == m)
							continue;
						inds[iter++] = (struct index){m,m,p,q};
					}
				}
			}

			for (u32 A = 0; A < res.component_count; ++A) {
//#pragma omp parallel for reduction(+: E3)
				for (u32 i = 0; i < len; ++i) {
					//sbmf_log_info("\t(%u,%u),(%u,%u)", inds[i].m,  inds[i].n, inds[i].p, inds[i].q);
					f64 v_AA_mmpq = V(&pt, A,A, inds[i].m, inds[i].n, inds[i].p, inds[i].q);

					/* Two double substitutions */
					// f64 me0 = (1.0/sqrt(2.0)) * G0(&pt,A,A) * v_AA_mmpq;
					f64 me0 = sqrt(2.0) * G0(&pt,A,A) * v_AA_mmpq;

					/* One double substitution */
					f64 me1 	= rs_2nd_order_me(&pt, A,A, inds[i].m, inds[i].n);
					f64 Ediff0  = rs_2nd_order_ediff(&pt, A,A, inds[i].m, inds[i].n);

					f64 me2 	= rs_2nd_order_me(&pt, A,A, inds[i].p, inds[i].q);
					f64 Ediff1  = rs_2nd_order_ediff(&pt, A,A, inds[i].p, inds[i].q);

					E3 += (me1*me0*me2)/(Ediff0*Ediff1);
				}
			}
		}
		sbmf_log_info("\tE3: %lf -- contrib: %lf", E3, E3-E3_before);

		/* mp,mq */
		/* states_to_include*(states_to_include - 1)*(states_to_include-2)/2 */
		sbmf_log_info("mp,mq");
		E3_before = E3;
		{
			const u32 len = (states_to_include-1)*(states_to_include-2)*(states_to_include-3)/2;
			u32 inds[len][4];
			u32 iter = 0;
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 p = 1; p < states_to_include; ++p) {
					if (p == m)
						continue;
					for (u32 q = p+1; q < states_to_include; ++q) {
						if (q == m)
							continue;
						inds[iter][0] = m;
						inds[iter][1] = p;
						inds[iter][2] = m;
						inds[iter][3] = q;
						iter++;
					}
				}
			}

			for (u32 A = 0; A < res.component_count; ++A) {
//#pragma omp parallel for reduction(+: E3)
				for (u32 i = 0; i < len; ++i) {
					//sbmf_log_info("\t(%u,%u),(%u,%u)", inds[i][0],  inds[i][1], inds[i][2], inds[i][3]);
					/* p0q0 */
					f64 v_AA_p0q0 = V(&pt, A,A, inds[i][1],0,inds[i][3],0);
					/* mpmq */
					f64 v_AA_mpmq = V(&pt, A,A, inds[i][0], inds[i][1], inds[i][2], inds[i][3]);

					/* Two double substitutions */
					f64 me0 = G0(&pt,A,A) * ((particle_count[A]-3)*v_AA_p0q0 + 2*v_AA_mpmq);

					/* One double substitution */
					f64 me1 	= rs_2nd_order_me(&pt,    A,A, inds[i][0], inds[i][1]);
					f64 Ediff0  = rs_2nd_order_ediff(&pt, A,A, inds[i][0], inds[i][1]);
					f64 me2 	= rs_2nd_order_me(&pt,    A,A, inds[i][2], inds[i][3]);
					f64 Ediff1  = rs_2nd_order_ediff(&pt, A,A, inds[i][2], inds[i][3]);
					E3 += (me1*me0*me2)/(Ediff0*Ediff1);
				}
			}
		}
		sbmf_log_info("\tE3: %lf -- contrib: %lf", E3, E3-E3_before);

		/* mn,pq */
		/* states_to_include choose 4 */
		sbmf_log_info("mn,pq");
		E3_before = E3;
		{
			const u32 len = n_choose_k(states_to_include-1, 4);
			u32 inds[len][4];
			u32 iter = 0;
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 n = m+1; n < states_to_include; ++n) {
					for (u32 p = n+1; p < states_to_include; ++p) {
						for (u32 q = p+1; q < states_to_include; ++q) {
							inds[iter][0] = m;
							inds[iter][1] = n;
							inds[iter][2] = p;
							inds[iter][3] = q;
							iter++;
						}
					}
				}
			}

			for (u32 A = 0; A < res.component_count; ++A) {
//#pragma omp parallel for reduction(+: E3)
				for (u32 i = 0; i < len; ++i) {
					//sbmf_log_info("\t(%u,%u),(%u,%u)", inds[i][0],  inds[i][1], inds[i][2], inds[i][3]);
					f64 v_AA_mnpq = V(&pt, A,A, inds[i][0], inds[i][1], inds[i][2], inds[i][3]);

					/* Two double substitutions */
					f64 me0 = 2*G0(&pt,A,A)*v_AA_mnpq;

					/* One double substitution */
					f64 me1 	= rs_2nd_order_me(&pt,    A,A, inds[i][0], inds[i][1]);
					f64 Ediff0  = rs_2nd_order_ediff(&pt, A,A, inds[i][0], inds[i][1]);
					f64 me2 	= rs_2nd_order_me(&pt,    A,A, inds[i][2], inds[i][3]);
					f64 Ediff1  = rs_2nd_order_ediff(&pt, A,A, inds[i][2], inds[i][3]);

					E3 += (me1*me0*me2)/(Ediff0*Ediff1);
				}
			}
		}
		sbmf_log_info("\tE3: %lf -- contrib: %lf", E3, E3-E3_before);

#if 1
		/* AmBn,AmBn */
		sbmf_log_info("AmBn,AmBn; AmBm,AmBm");
		sbmf_log_info("E3 before: %e", E3);
		for (u32 A = 0; A < res.component_count; ++A) {
			for (u32 B = A+1; B < res.component_count; ++B) {
//#pragma omp parallel for reduction(+: E3)
				for (u32 m = 1; m < states_to_include; ++m) {
					for (u32 n = 1; n < states_to_include; ++n) {
						sbmf_log_info("\t(%u,%u) -- (%u,%u),(%u,%u)", A,B, m,n,m,n);
						f64 v_AA_0000 = V(&pt, A,A, 0,0,0,0);
						f64 v_BB_0000 = V(&pt, B,B, 0,0,0,0);
						f64 v_AA_m0m0 = V(&pt, A,A, m,0,m,0);
						f64 v_BB_n0n0 = V(&pt, B,B, n,0,n,0);

						f64 v_AB_0000 = V(&pt, A,B, 0,0,0,0);
						f64 v_AB_mnmn = V(&pt, A,B, m,n,m,n);
						f64 v_AB_m0m0 = V(&pt, A,B, m,0,m,0);
						f64 v_AB_0n0n = V(&pt, A,B, 0,n,0,n);

						f64 me0 = 0;
						me0 += 0.5*G0(&pt,A,A)*(-particle_count[A]*(particle_count[A]-1)*v_AA_0000
								+ 2*(particle_count[A]-1)*v_AA_m0m0);
						me0 += 0.5*G0(&pt,B,B)*(-particle_count[B]*(particle_count[B]-1)*v_BB_0000
								+ 2*(particle_count[B]-1)*v_BB_n0n0);
						me0 += G0(&pt,A,B)*(
								(1-particle_count[A]*particle_count[B])*v_AB_0000
								-v_AB_0n0n
								-v_AB_m0m0
								+v_AB_mnmn
								);


						/* One double substituion */
						f64 me1 	= rs_2nd_order_me(&pt,    A,B, m,n);
						f64 Ediff0  = rs_2nd_order_ediff(&pt, A,B, m,n);

						E3 += (me1*me0*me1)/(Ediff0*Ediff0);
					}
				}

			}
		}

		/* AmBn,AmBp */
		sbmf_log_info("AmBn,AmBp");
		sbmf_log_info("E3 before: %e", E3);
		for (u32 A = 0; A < res.component_count; ++A) {
			for (u32 B = 0; B < res.component_count; ++B) {
				if (B == A)
					continue;

//#pragma omp parallel for reduction(+: E3)
				for (u32 m = 1; m < states_to_include; ++m) {
					for (u32 n = 1; n < states_to_include; ++n) {
						if (n == m) continue;
						for (u32 p = n+1; p < states_to_include; ++p) {
							if (p == m) continue;
							sbmf_log_info("\t(%u,%u) -- (%u,%u),(%u,%u)", A,B, m,n,m,p);
							f64 v_BB_0n0p = V(&pt, B,B, 0,n,0,p);
							f64 v_AB_mnmp = V(&pt, A,B, m,n,m,p);
							f64 v_AB_0n0p = V(&pt, A,B, 0,n,0,p);

							f64 me0 = G0(&pt,B,B)*(particle_count[B]-1)*v_BB_0n0p
								+ G0(&pt,A,B)*(v_AB_mnmp - v_AB_0n0p);

							/* One double substition */
							f64 me1 	= rs_2nd_order_me(&pt,    A,B, m,n);
							f64 Ediff0  = rs_2nd_order_ediff(&pt, A,B, m,n);
							f64 me2 	= rs_2nd_order_me(&pt,    A,B, m,p);
							f64 Ediff1  = rs_2nd_order_ediff(&pt, A,B, m,p);

							E3 += (me1*me0*me2)/(Ediff0*Ediff1);
						}
					}
				}
			}
		}


		/* AmBn,ApBq */
		sbmf_log_info("AmBn,ApBq");
		sbmf_log_info("E3 before: %e", E3);
		for (u32 A = 0; A < res.component_count; ++A) {
			for (u32 B = A+1; B < res.component_count; ++B) {

		/* states_to_include * (states_to_include - 1) */
//#pragma omp parallel for reduction(+: E3)
				for (u32 m = 1; m < states_to_include; ++m) {
					for (u32 n = 1; n < states_to_include; ++n){
						if (n == m) continue;
						for (u32 p = 1; p < states_to_include; ++p) {
							if (p == n || p == m) continue;
							for (u32 q = 1; q < states_to_include; ++q) {
								if (q == p || q == m || q == n) continue;
								sbmf_log_info("\t(%u,%u) -- (%u,%u),(%u,%u)", A,B, m,n,p,q);

								/* Making sure we avoid mnpq <-> pqmn */
								if (m*10+n > p*10+q) continue;

								f64 v_AB_mnpq = V(&pt, A,B, m,n,p,q);

								f64 me0 = G0(&pt, A,B) * v_AB_mnpq;

								/* One double substitution */
								f64 me1 	= rs_2nd_order_me(&pt,    A,B, m,n);
								f64 Ediff0  = rs_2nd_order_ediff(&pt, A,B, m,n);
								f64 me2 	= rs_2nd_order_me(&pt,    A,B, p,q);
								f64 Ediff1  = rs_2nd_order_ediff(&pt, A,B, p,q);

								E3 += (me1*me0*me2)/(Ediff0*Ediff1);
							}
						}
					}
				}

			}
		}
#endif

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
