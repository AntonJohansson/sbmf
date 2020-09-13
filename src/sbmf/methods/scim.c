#include "find_groundstate.h"
#include <sbmf/sbmf.h>
#include <sbmf/math/find_eigenpairs.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/methods/quadgk_vec.h>
#include <sbmf/debug/profile.h>

#include <assert.h> /* Not the correct way to handle this */

struct integrand_params {
	u32 n[2];
	c64* coeffs;
	u32 coeff_count;
	gss_potential_vec_func* pot;
};

struct guess_integrand_params {
	u32 n;
	gss_guess_vec_func* guess;
};

static void scim_integrand(f64* out, f64* in, u32 len, void* data) {
	struct integrand_params* params = data;

	f64 eig1[len];
	f64 eig2[len];
	ho_eigenfunction_vec(params->n[0], eig1, in, len);
	ho_eigenfunction_vec(params->n[1], eig2, in, len);

	c64 sample[len];
	hob_sample_vec(params->coeffs, params->coeff_count, sample, in, len);

	f64 pot[len];
	params->pot(pot, in, sample, len);

	for (u32 i = 0; i < len; ++i) {
		out[i] = eig1[i]*eig2[i]*pot[i];
	}

	/*
	if (isnan(res)) {
		log_error("res is nan!");
		log_error("eig1: %lf", eig1);
		log_error("eig2: %lf", eig2);
		log_error("pot: %lf", pot);
		log_error("x: %lf", x);
	}
	*/
}

static void scim_guess_integrand(f64* out, f64* in, u32 len, void* data) {
	struct guess_integrand_params* params = data;

	f64 eig[len];
	ho_eigenfunction_vec(params->n, eig, in, len);

	c64 guess[len];
	params->guess(guess, in, len);

	for (u32 i = 0; i < len; ++i) {
		out[i] = eig[i]*eig[i]*guess[i];
	}
}

struct gss_result scim(struct scim_settings settings, gss_potential_vec_func* potential, gss_guess_vec_func* guess) {
	const u32 N = settings.num_basis_functions;

	struct gss_result res = {
		.wavefunction = (c64*) sbmf_stack_push(N*sizeof(c64)),
	};
	c64* old_wavefunction = (c64*) sbmf_stack_push(N*sizeof(c64));

	integration_settings int_settings = {
		.gk = gk7,
		.abs_error_tol = 1e-10,
		.rel_error_tol = 1e-10,
		.max_evals = 1e4,
	};

	/* Setup the initial guess */
	/*res.wavefunction[0] = 1;
	for (u32 i = 1; i < N; ++i) {
		res.wavefunction[i] = 0;
	}*/
	struct guess_integrand_params guess_params = {
		.guess = guess,
	};
	int_settings.userdata = &guess_params;
	for (u32 i = 0; i < N; ++i) {
		guess_params.n = i;
		integration_result ires = quadgk_vec(scim_guess_integrand, -INFINITY, INFINITY, int_settings);
		if (!ires.converged)
			log_error("integration failed for %d", i);
		assert(ires.converged);
		res.wavefunction[i] = ires.integral;
	}

	hermitian_bandmat T = construct_ho_kinetic_matrix(N);
	hermitian_bandmat H = T;
	H.base.data = (c64*) sbmf_stack_push(N*N*sizeof(c64));

	struct integrand_params params = {
		.pot = potential,
		.coeff_count = N,
		.coeffs = res.wavefunction,
	};

	int_settings.userdata = &params;


	for (; res.iterations < settings.max_iterations; ++res.iterations) {
		memcpy(old_wavefunction, res.wavefunction, N*sizeof(c64));

		log_info("Staring iteration: %d -- error: %lf", res.iterations, res.error);

		/* Construct standard hamiltonian */
		{
			PROFILE_BEGIN("Constructing H");
			for (u32 r = 0; r < T.size; ++r) {
				for (u32 c = r; c < T.size; ++c) {
					params.n[0] = r;
					params.n[1] = c;
					integration_result res = quadgk_vec(scim_integrand, -INFINITY, INFINITY, int_settings);
					if (!res.converged)
						log_error("integration failed for %d,%d", r,c);
					assert(res.converged);

					u32 i = (H.size-1)*(H.size-(c-r)) + r;
					H.base.data[i] = T.base.data[i] + res.integral;
				}
			}
			PROFILE_END("Constructing H");

			assert(mat_is_valid(H.base));
		}

		/* Solve for first eigenvector (ground state) */
		PROFILE_BEGIN("Eigenproblem solving");
		struct eigen_result eres = find_eigenpairs_sparse(H, 1, EV_SMALLEST_RE);
		PROFILE_END("Eigenproblem solving");

		/* Normalize and copy to result */
		{
			PROFILE_BEGIN("Normalize+copy");
			f64 sum = 0.0;
			for (u32 i = 0; i < N; ++i) {
				res.wavefunction[i] = eres.eigenvectors[i];
				f64 tmp = cabs(eres.eigenvectors[i]);
				sum += tmp*tmp;
			}

			f64 scaling = 1.0/sqrt(sum);
			for (u32 i = 0; i < N; ++i) {
				res.wavefunction[i] *= scaling;
			}
			PROFILE_END("Normalize+copy");
		}

		/* Calculate error */
		{
			PROFILE_BEGIN("Calculate error");
			f64 sum = 0.0;
			for (u32 i = 0; i < N; ++i) {
				f64 diff = cabs(res.wavefunction[i]) - cabs(old_wavefunction[i]);
				sum += diff*diff;
			}
			res.error = sqrt(sum);

			PROFILE_END("Calculate error");
			if (res.error < settings.error_tol)
				break;
		}

		/* Call debug callback if requested by user */
		if (settings.measure_every > 0 && res.iterations % settings.measure_every == 0) {
			if (settings.dbgcallback) {
				/*
				 * Note: old_wavefunction has already been used at this point, we can
				 * thus use it as a temporary buffer.
				 */
				settings.dbgcallback(settings, res.wavefunction);
			}
		}
	}

	return res;
}
