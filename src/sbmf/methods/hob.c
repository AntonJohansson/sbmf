#include "find_groundstate.h"
#include <sbmf/sbmf.h>
#include <sbmf/math/find_eigenpairs.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/methods/quadgk.h>
#include <sbmf/debug/profile.h>

#include <assert.h> /* Not the correct way to handle this */

static inline f64 hob_integrand(f64 x, void* data) {
	PROFILE_BEGIN("integrand -- total");
	hob_integrand_params* params = data;

	PROFILE_BEGIN("integrand -- eigs");
	f64 eig1 = ho_eigenfunction((i32[]){params->n[0]}, &x, 1);
	f64 eig2 = ho_eigenfunction((i32[]){params->n[1]}, &x, 1);
	PROFILE_END("integrand -- eigs");

	//PROFILE_BEGIN("integrand sample");
	//f64 sample = hob_sample(params->coeffs, params->coeff_count, x);
	//PROFILE_END("integrand sample");
	f64 sample = 0;

	PROFILE_BEGIN("integrand -- pot");
	f64 pot = params->pot(&x,1, sample);
	PROFILE_END("integrand -- pot");

	f64 res = (eig1*eig2)*pot;
	PROFILE_END("integrand -- total");

	if (isnan(res)) {
		log_info("eig1: %lf", eig1);
		log_info("eig2: %lf", eig2);
		log_info("pot: %lf", pot);
		log_info("x: %lf", x);
	}
	return res;
}

gss_result hob(gss_settings settings, gss_potential_func* potential, gss_guess_func* guess) {
	const u32 N = 128;

	gss_result res = {
		.settings = settings,
		.wavefunction = (c64*) sbmf_stack_push(N*sizeof(c64)),
	};
	c64* old_wavefunction = (c64*) sbmf_stack_push(N*sizeof(c64));

	res.wavefunction[0] = 1;
	for (u32 i = 1; i < N; ++i) {
		res.wavefunction[i] = 0;
	}

	hermitian_bandmat T = construct_ho_kinetic_matrix(N);
	hermitian_bandmat H = T;
	H.base.data = (c64*) sbmf_stack_push(N*N*sizeof(c64));

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
					if (!res.converged)
						log_error("integration failed for %d,%d", r,c);
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
		struct eigen_result eres = find_eigenpairs_sparse(H, 1, EV_SMALLEST_RE);
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

/*
#define FACTORIAL_MAX_N 20

static u64 factorial(u32 n) {
	assert(n <= FACTORIAL_MAX_N); // largest factorial supported by u64

	u64 res = 1;
	while (n > 1)
		res *= n--;

	return res;
}

static inline f64 ho_eigenfunction(i32 states[], f64 point[], i32 dims) {
	f64 prod = 1.0;
	static const f64 pi_factor = 1.0/pow(M_PI,0.25);
	for (i32 i = 0; i < dims; ++i) {
		i32 n = states[i];
		f64 x = point[i];

		f64 factorial_value = (n <= FACTORIAL_MAX_N) ? factorial(n) : sqrt(M_2_PI*n) * pow(n/exp(1), n);
		f64 normalization_factor = exp(-x*x/2.0) * pi_factor / sqrt(pow(2,n) * factorial_value);

		{
			f64 H0 = normalization_factor*1;
			f64 H1 = normalization_factor*2*x;
			f64 HN = 0.0;

			if (n == 0) {
				HN = H0;
			} else if (n == 1) {
				HN = H1;
			} else {
				for (i32 i = 2; i <= n; ++i) {
					HN = 2*x*H1 - 2*(i-1)*H0;
					H0 = H1;
					H1 = HN;
				}
			}

			prod *= HN;
		}
	}
	return prod;
}
*/
