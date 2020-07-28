#include "groundstate_solver.h"
#include <sbmf/common/common.h>
#include <sbmf/common/eigenproblem.h>
#include <sbmf/basis/harmonic_oscillator.h>
#include <sbmf/numerical_integration/quadgk.h>
#include <sbmf/common/profile.h>
#include <stdlib.h>
#include <assert.h>

static inline f64 hob_integrand(f64 x, void* data) {
	PROFILE_BEGIN("integrand");
	hob_integrand_params* params = data;

	PROFILE_BEGIN("integrand eigs");
	f64 eig1 = ho_eigenfunction((i32[]){params->n[0]}, &x, 1);
	f64 eig2 = ho_eigenfunction((i32[]){params->n[1]}, &x, 1);
	PROFILE_END("integrand eigs");

	//PROFILE_BEGIN("integrand sample");
	//f64 sample = hob_sample(params->coeffs, params->coeff_count, x);
	//PROFILE_END("integrand sample");
	f64 sample = 0;

	PROFILE_BEGIN("integrand pot");
	f64 pot = params->pot(&x,1, sample);
	PROFILE_END("integrand pot");

	f64 res = (eig1*eig2)*pot;
	PROFILE_END("integrand");
	return res;
}

gss_result hob(gss_settings settings, gss_potential_func* potential, gss_guess_func* guess) {
	const u32 N = 64;

	gss_result res = {
		.settings = settings,
		.wavefunction = (c64*) sa_push(_sbmf.main_stack, N*sizeof(c64)),
	};
	c64* old_wavefunction = (c64*)sa_push(_sbmf.main_stack, N*sizeof(c64));

	res.wavefunction[0] = 1;
	for (u32 i = 1; i < N; ++i) {
		res.wavefunction[i] = 0;
	}

	hermitian_bandmat T = construct_ho_kinetic_matrix(N);
	hermitian_bandmat H = T;
	H.base.data = (c64*) sa_push(_sbmf.main_stack, N*N*sizeof(c64));

	hob_integrand_params params = {
		.pot = potential,
		.coeff_count = N,
		.coeffs = res.wavefunction,
	};
	integration_settings int_settings = {
		.gk = gk7,
		.abs_error_tol = 1e-10,
		.rel_error_tol = 1e-10,
		.max_evals = 1e4,
		.userdata = &params
	};

	for (; res.iterations < settings.max_iterations; ++res.iterations) {
		log_info("%d", res.iterations);
		memcpy(old_wavefunction, res.wavefunction, N*sizeof(c64));
		// Construct standars hamiltonian
		{
			PROFILE_BEGIN("H");
			for (u32 r = 0; r < T.size; ++r) {
				for (u32 c = r; c < T.size; ++c) {
					PROFILE_BEGIN("H iter");
					params.n[0] = r;
					params.n[1] = c;
					integration_result res = quadgk(hob_integrand, -INFINITY, INFINITY, int_settings);
					assert(res.converged);

					u32 i = (H.size-1)*(H.size-(c-r)) + r;
					H.base.data[i] = T.base.data[i] + res.integral;
					PROFILE_END("H iter");
				}
			}
			PROFILE_END("H");

			assert(mat_is_valid(H.base));
		}
		// Solve for first eigenvector (ground state)
		eig_result eres = eig_sparse_bandmat(H, 1, EV_SMALLEST_RE);
		// Normalize and copy to result
		{
			f64 sum = 0.0;
			for (u32 i = 0; i < N; ++i) {
				res.wavefunction[i] = eres.eigenvectors[i];
				f64 tmp = cabs(eres.eigenvectors[i]);
				sum += tmp*tmp;
			}

			//f64 scaling = 1.0/sqrt(sum);
			//for (u32 i = 0; i < N; ++i) {
			//	res.wavefunction[i] *= scaling;
			//}
		}
		// Calculate error
		{
			f64 sum = 0.0;
			for (u32 i = 0; i < N; ++i) {
				f64 diff = cabs(res.wavefunction[i]) - cabs(old_wavefunction[i]);
				sum += diff*diff;
			}
			res.error = sqrt(sum);

			log_info("error: %e", res.error);
			if (res.error < settings.error_tol)
				break;
		}
	}

	return res;
}
