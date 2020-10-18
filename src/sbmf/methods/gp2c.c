#include "find_groundstate.h"
#include <sbmf/sbmf.h>
#include <sbmf/math/find_eigenpairs.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/methods/quadgk_vec.h>
#include <sbmf/math/functions.h>
#include <assert.h> /* Not the correct way to handle this */
#include <omp.h>

struct integrand_params {
	u32 n[2];
	c64* coeff_a;
	c64* coeff_b;
	u32 coeff_count;
	gp2c_operator_func* op;
};

static void log_integration_result(integration_result res) {
	log_info("integral: %e", res.integral);
	log_info("error: %e", res.error);
	log_info("performed evals: %d", res.performed_evals);
	log_info("converged: %s", (res.converged) ? "yes" : "no");
}

static void scim_integrand(f64* out, f64* in, u32 len, void* data) {
	struct integrand_params* params = data;

	f64 eig1[len];
	f64 eig2[len];
	ho_eigenfunction_vec(params->n[0], eig1, in, len);
	ho_eigenfunction_vec(params->n[1], eig2, in, len);

	c64 sample_a[len];
	hob_sample_vec(params->coeff_a, params->coeff_count, sample_a, in, len);

	c64 sample_b[len];
	hob_sample_vec(params->coeff_b, params->coeff_count, sample_b, in, len);

	f64 op[len];
	params->op(op, in, sample_a, sample_b, len);

	for (u32 i = 0; i < len; ++i) {
		out[i] = eig1[i]*eig2[i]*op[i];
	}

#if 0
	for (u32 i = 0; i < len; ++i) {
		if (isnan(out[i])) {
			log_error("out[%d] is nan!", i);
			log_error("eig1: %lf", eig1[i]);
			log_error("eig2: %lf", eig2[i]);
			log_error("pot: %lf", pot[i]);
			log_error("in: %lf", in[i]);
		}
	}
#endif
}

struct gp2c_result gp2c(struct gp2c_settings settings, gp2c_operator_func* op_a, gp2c_operator_func* op_b) {
	const u32 N = settings.num_basis_functions;

	struct gp2c_result res = {
		.coeff_a = (c64*) sbmf_stack_push(N*sizeof(c64)),
		.coeff_b = (c64*) sbmf_stack_push(N*sizeof(c64)),
		.hamiltonian_a = complex_hermitian_bandmat_new_zero(N,N),
		.hamiltonian_b = complex_hermitian_bandmat_new_zero(N,N),
	};
	c64* old_coeff_a = (c64*) sbmf_stack_push(N*sizeof(c64));
	c64* old_coeff_b = (c64*) sbmf_stack_push(N*sizeof(c64));

	integration_settings int_settings = {
		.gk = gk7,
		.abs_error_tol = 1e-10,
		.rel_error_tol = 1e-10,
		.max_evals = settings.max_iterations,
	};

	/* Set initial guess to the ho groundstate */
	{
		if (settings.guess_a) {
			settings.guess_a(res.coeff_a, N);
		} else {
			for (u32 i = 0; i < N; ++i) {
				res.coeff_a[i] = 0;
			}
			res.coeff_a[0] = 1;
		}
		if (settings.guess_b) {
			settings.guess_b(res.coeff_b, N);
		} else {
			for (u32 i = 0; i < N; ++i) {
				res.coeff_b[i] = 0;
			}
			res.coeff_b[0] = 1;
		}
	}

	struct integrand_params params_a = {
		.coeff_count = N,
		.coeff_a = res.coeff_a,
		.coeff_b = res.coeff_b,
		.op = op_a,
	};

	struct integrand_params params_b = {
		.coeff_count = N,
		.coeff_a = res.coeff_a,
		.coeff_b = res.coeff_b,
		.op = op_b,
	};

	/* Precompute indices for matrix elements
	 * 	This way we get a single array looping
	 * 	over the matrix which is easier to
	 * 	parallelize well.
	 * */
	u32 matrix_element_rows[N*(N+1)/2];
	u32 matrix_element_cols[N*(N+1)/2];
	{
		u32 matrix_element_index = 0;
		COMPLEX_HERMITIAN_BANDMAT_FOREACH(res.hamiltonian_a, r,c) {
			matrix_element_rows[matrix_element_index] = r;
			matrix_element_cols[matrix_element_index] = c;
			matrix_element_index++;
		}
	}

	/* Do the actual iterations */
	for (; res.iterations < settings.max_iterations; ++res.iterations) {
		memcpy(old_coeff_a, res.coeff_a, N*sizeof(c64));
		memcpy(old_coeff_b, res.coeff_b, N*sizeof(c64));

		/* Construct standard hamiltonian for a */
		{
#pragma omp parallel for firstprivate(params_a, int_settings) shared(res)
			for (u32 i = 0; i < N*(N+1)/2; ++i) {
				u32 r = matrix_element_rows[i];
				u32 c = matrix_element_cols[i];
				int_settings.userdata = &params_a;
				params_a.n[0] = r;
				params_a.n[1] = c;

				integration_result int_res = quadgk_vec(scim_integrand, -INFINITY, INFINITY, int_settings);

				if (!int_res.converged) {
					log_error("integration failed for %d,%d", r,c);
					log_integration_result(int_res);
				}
				assert(int_res.converged);

				u32 i = complex_hermitian_bandmat_index(res.hamiltonian_a, r,c);
				res.hamiltonian_a.data[i] = int_res.integral;
				if (r == c)
					res.hamiltonian_a.data[i] += ho_eigenvalue((i32[]){r},1);
			}

			assert(complex_hermitian_bandmat_is_valid(res.hamiltonian_a));
		}

		/* Construct standard hamiltonian for b */
		{
#pragma omp parallel for firstprivate(params_b, int_settings) shared(res)
			for (u32 i = 0; i < N*(N+1)/2; ++i) {
				u32 r = matrix_element_rows[i];
				u32 c = matrix_element_cols[i];
				int_settings.userdata = &params_b;
				params_b.n[0] = r;
				params_b.n[1] = c;

				integration_result int_res = quadgk_vec(scim_integrand, -INFINITY, INFINITY, int_settings);

				if (!int_res.converged) {
					log_error("integration failed for %d,%d", r,c);
					log_integration_result(int_res);
				}
				assert(int_res.converged);

				u32 i = complex_hermitian_bandmat_index(res.hamiltonian_b, r,c);
				res.hamiltonian_b.data[i] = int_res.integral;
				if (r == c)
					res.hamiltonian_b.data[i] += ho_eigenvalue((i32[]){r},1);
			}

			assert(complex_hermitian_bandmat_is_valid(res.hamiltonian_b));
		}

		/* Solve for first eigenvector (ground state) */
		struct eigen_result eres_a = find_eigenpairs_sparse(res.hamiltonian_a, 1, EV_SMALLEST_RE);
		struct eigen_result eres_b = find_eigenpairs_sparse(res.hamiltonian_b, 1, EV_SMALLEST_RE);

		/* Copy energies */
		res.energy_a = eres_a.eigenvalues[0];
		res.energy_b = eres_b.eigenvalues[0];

		/* Normalize and copy to result */
		{
			c64_normalize(res.coeff_a, eres_a.eigenvectors, N);
			c64_normalize(res.coeff_b, eres_b.eigenvectors, N);
		}

		/* Callback */
		if (settings.post_normalize_callback) {
			settings.post_normalize_callback(res.coeff_a, res.coeff_b, N);
		}

		/* Calculate error */
		{
			f64 sum_a = 0.0;
#pragma omp parallel for shared(res, old_coeff_a) reduction(+: sum_a)
			for (u32 i = 0; i < N; ++i) {
				f64 diff_a = cabs(res.coeff_a[i]) - cabs(old_coeff_a[i]);
				sum_a += diff_a*diff_a;
			}
			res.error_a = sqrt(sum_a);

			f64 sum_b = 0.0;
#pragma omp parallel for shared(res, old_coeff_b) reduction(+: sum_b)
			for (u32 i = 0; i < N; ++i) {
				f64 diff_b = cabs(res.coeff_b[i]) - cabs(old_coeff_b[i]);
				sum_b += diff_b*diff_b;
			}
			res.error_b = sqrt(sum_b);
		}

		/* Break condition */
		log_info("Finished iteration: %d with energyies: a: %.2e, b: %.2e, ana error of a: %.2e, b: %.2e", res.iterations, res.energy_a, res.energy_b, res.error_a, res.error_b);
		if (res.error_a < settings.error_tol && res.error_b < settings.error_tol)
			break;

#if 0
		/* Call debug callback if requested by user */
		if (settings.measure_every > 0 && res.iterations % settings.measure_every == 0) {
			if (settings.dbgcallback) {
				settings.dbgcallback(settings, res.wavefunction);
			}
		}
#endif
	}

	return res;
}
