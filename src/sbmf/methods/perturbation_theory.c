#include <sbmf/methods/perturbation_theory.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/math/find_eigenpairs.h>
#include <sbmf/math/functions.h>
#include <sbmf/sbmf.h>
#include <assert.h>

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

f64 V(const u32 coeff_count, f64 i[static coeff_count],
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

	struct integration_result res = quadgk_vec(V_integrand, -INFINITY, INFINITY, settings);
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
					E2 += me*me/(2*ENERGY(A,0) - (ENERGY(A,i) + ENERGY(A,j)));
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
						f64 me = G0(A,B) * sqrt(particle_count[i] * particle_count[j])
								* V(L, &PHI(A,i), &PHI(B,j), &PHI(A,0), &PHI(B,0));
						E2 += me*me/(ENERGY(A,0) + ENERGY(B,0) - (ENERGY(A,i) + ENERGY(B,j)));
					}
				}
			}
		}
	}

	return (struct pt_result) {
		.E0 = E0,
		.E1 = E1,
		.E2 = E2,
	};
}
