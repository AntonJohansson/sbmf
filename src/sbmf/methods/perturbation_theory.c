/*
 * Holds all information needed to do the PT,
 * easy to pass around. Basicly params to
 * pt_rayleigh_schroedinger
 */
struct pt_settings {
	struct nlse_result* res;
	f64* g0;
	i32* particle_count;

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

static f64 V(struct pt_settings* pt, u32 A, u32 B, u32 i, u32 j, u32 k, u32 l) {
	/* Check if we have computed this integral already */
	struct V_params p = {
		.coeff_count = pt->L,
		.i = PHI(pt, A, i),
		.j = PHI(pt, B, j),
		.k = PHI(pt, A, k),
		.l = PHI(pt, B, l)
	};

    struct integration_settings settings = {
        .abs_error_tol = 1e-10,
        .rel_error_tol = 1e-7,
		.gk = gk15,
        .max_evals = 1e5,
		.userdata = &p,
    };

	struct integration_result res = quadgk(V_integrand, -INFINITY, INFINITY, settings);
	assert(res.converged);

	return res.integral;
}

/*
 * Helper functions since these will be calculated a lot
 */

static f64 rs_2nd_order_me(struct pt_settings* pt, u32 A, u32 B, u32 i, u32 j) {
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

static f64 rs_2nd_order_ediff(struct pt_settings* pt, u32 A, u32 B, u32 i, u32 j) {
	return E(pt,A,0) + E(pt,B,0) - E(pt,A,i) - E(pt,B,j);
}

/*
 * Main function for Rayleigh-Schrodinger perturbation theory
 */

struct pt_result rayleigh_schroedinger_pt(struct nlse_result res, f64* g0, i32* particle_count) {
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

	struct pt_settings pt = {
		.res = &res,
		.g0 = g0,
		.particle_count = particle_count,
		.states = states,
		.N = states_to_include,
		.L = res.coeff_count,
	};

	/*
	 * Macros are not local, but they are only used here...
	 * PHI returns the coeffs. of state 'state' in component
	 * 'component'. ENERGY returns the energy of that same
	 * state. G0 return the interaction strength of component
	 * A with respect to component B, the g0 array is assumed
	 * to be symmetric.
	 */
#define ENERGY(component, state) states[component].eigenvalues[state]

	/* zeroth order PT, computes <0|H0|0> */
	f64 E0 = 0.0;
	{
		for (u32 A = 0; A < res.component_count; ++A) {
			E0 += particle_count[A] * E(&pt, A, 0);
		}
	}

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
						f64 me = rs_2nd_order_me(&pt, A,A, i,j);
						f64 Ediff = rs_2nd_order_ediff(&pt, A,B, i,j);
						E2 += me*me/(Ediff);
						E3_last_term += me*me/(Ediff*Ediff);
					}
				}
			}
		}
	}

	f64 E3 = 0.0;
#if 1
	{
		/*
		 * Loop over all <AmAn|V|ApAq>,
		 * requires all unique combinations (m,n,p,q)
		 */

		/* mm,mm */
		sbmf_log_info("mm,mm");
		for (u32 A = 0; A < res.component_count; ++A) {
#pragma omp parallel for reduction(+: E3)
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
				//for (u32 B = 0; B < res.component_count; ++B) {
				//	if (B == A)
				//		continue;

				//	f64 v_BA_0000 = V(&pt, B,A, 0,0,0,0);
				//	me0 -= G0(&pt,A,B)*particle_count[A]*particle_count[B]*v_BA_0000;
				//}

				/* One double substitution */
				f64 me1 = rs_2nd_order_me(&pt, A,A, m,m);
				f64 Ediff = rs_2nd_order_ediff(&pt, A,A, m,m);

				E3 += (me1*me0*me1)/(Ediff*Ediff);
			}
		}

		/* mn,mn */
		sbmf_log_info("mn,mn");
		for (u32 A = 0; A < res.component_count; ++A) {
#pragma omp parallel for reduction(+: E3)
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
					//for (u32 B = 0; B < res.component_count; ++B) {
					//	if (B == A)
					//		continue;

					//	f64 v_BA_0000 = V(&pt, B,A, 0,0,0,0);
					//	me0 -= G0(&pt,A,B)*particle_count[A]*particle_count[B]*v_BA_0000;
					//}

					/* One double substitution */
					f64 me1 = rs_2nd_order_me(&pt, A,A, m,n);
					f64 Ediff = rs_2nd_order_ediff(&pt, A,A, m,n);
					E3 += (me1*me0*me1)/(Ediff*Ediff);
				}
			}
		}

		/* mm,nn */
		sbmf_log_info("mm,nn");
		for (u32 A = 0; A < res.component_count; ++A) {
#pragma omp parallel for reduction(+: E3)
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

		/* mm,mp */
		sbmf_log_info("mm,mp");
		for (u32 A = 0; A < res.component_count; ++A) {
#pragma omp parallel for reduction(+: E3)
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

		/* mm,pq */
		sbmf_log_info("mm,pq");
		for (u32 A = 0; A < res.component_count; ++A) {
#pragma omp parallel for reduction(+: E3)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 p = 1; p < states_to_include; ++p) {
					if (p == m)
						continue;

					for (u32 q = p+1; q < states_to_include; ++q) {
						if (q == m)
							continue;
						//sbmf_log_info("\t(%u,%u),(%u,%u)", m,m,p,q);
						f64 v_AA_mmpq = V(&pt, A,A, m,m,p,q);

						/* Two double substitutions */
						//f64 me0 = (1.0/sqrt(2.0)) * G0(&pt,A,A) * v_AA_mmpq;
						f64 me0 = sqrt(2.0) * G0(&pt,A,A) * v_AA_mmpq;

						/* One double substitution */
						f64 me1 	= rs_2nd_order_me(&pt, A,A, m,m);
						f64 Ediff0  = rs_2nd_order_ediff(&pt, A,A, m,m);

						f64 me2 	= rs_2nd_order_me(&pt, A,A, p,q);
						f64 Ediff1  = rs_2nd_order_ediff(&pt, A,A, p,q);

						E3 += (me1*me0*me2)/(Ediff0*Ediff1);
					}
				}
			}
		}

		/* mp,mq */
		sbmf_log_info("mp,mq");
		for (u32 A = 0; A < res.component_count; ++A) {
#pragma omp parallel for reduction(+: E3)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 p = 1; p < states_to_include; ++p) {
					if (p == m)
						continue;

					for (u32 q = p+1; q < states_to_include; ++q) {
						if (q == m)
							continue;
						//sbmf_log_info("\t(%u,%u),(%u,%u)", m,p,m,q);
						f64 v_AA_p0q0 = V(&pt, A,A, p,0,q,0);
						f64 v_AA_mpmq = V(&pt, A,A, m,p,m,q);

						/* Two double substitutions */
						f64 me0 = G0(&pt,A,A) * ((particle_count[A]-3)*v_AA_p0q0 + 2*v_AA_mpmq);

						/* One double substitution */
						f64 me1 	= rs_2nd_order_me(&pt,    A,A, m,p);
						f64 Ediff0  = rs_2nd_order_ediff(&pt, A,A, m,p);
						f64 me2 	= rs_2nd_order_me(&pt,    A,A, m,q);
						f64 Ediff1  = rs_2nd_order_ediff(&pt, A,A, m,q);
						E3 += (me1*me0*me2)/(Ediff0*Ediff1);
					}
				}
			}
		}

		/* mn,pq */
		sbmf_log_info("mn,pq");
		for (u32 A = 0; A < res.component_count; ++A) {
#pragma omp parallel for reduction(+: E3)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 n = m+1; n < states_to_include; ++n) {
					for (u32 p = n+1; p < states_to_include; ++p) {
						for (u32 q = p+1; q < states_to_include; ++q) {
							//sbmf_log_info("\t(%u,%u),(%u,%u)", m,n,p,q);
							f64 v_AA_mnpq = V(&pt, A,A, m,n,p,q);

							/* Two double substitutions */
							f64 me0 = 2*G0(&pt,A,A)*v_AA_mnpq;

							/* One double substitution */
							f64 me1 	= rs_2nd_order_me(&pt,    A,A, m,n);
							f64 Ediff0  = rs_2nd_order_ediff(&pt, A,A, m,n);
							f64 me2 	= rs_2nd_order_me(&pt,    A,A, p,q);
							f64 Ediff1  = rs_2nd_order_ediff(&pt, A,A, p,q);

							E3 += (me1*me0*me2)/(Ediff0*Ediff1);
						}
					}
				}
			}
		}

#if 0
		/* AmBn,AmBn */
		sbmf_log_info("AmBn,AmBn; AmBm,AmBm");
		sbmf_log_info("E3 before: %e", E3);
		for (u32 A = 0; A < res.component_count; ++A) {
			for (u32 B = A+1; B < res.component_count; ++B) {
#pragma omp parallel for reduction(+: E3)
				for (u32 m = 1; m < states_to_include; ++m) {
					for (u32 n = 1; n < states_to_include; ++n) {
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

#pragma omp parallel for reduction(+: E3)
				for (u32 m = 1; m < states_to_include; ++m) {
					for (u32 n = 1; n < states_to_include; ++n) {
						if (n == m) continue;
						for (u32 p = n+1; p < states_to_include; ++p) {
							if (p == m) continue;
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

#pragma omp parallel for reduction(+: E3)
				for (u32 m = 1; m < states_to_include; ++m) {
					for (u32 n = 1; n < states_to_include; ++n){
						if (n == m) continue;
						for (u32 p = 1; p < states_to_include; ++p) {
							if (p == n || p == m) continue;
							for (u32 q = 1; q < states_to_include; ++q) {
								if (q == p || q == m || q == n) continue;

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
