#include <sbmf/methods/gp2c.h>
#include <sbmf/sbmf.h>
#include <sbmf/math/find_eigenpairs.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/math/functions.h>

#include <assert.h> /* Not the correct way to handle this */
#include <omp.h>

struct integrand_params {
	u32 n[2];
	u32 component_count;
	u32 coeff_count;
	f64* coeff;

	gp2c_operator_func* op;
	void* op_userdata;

	struct basis basis;
};

static void sbmf_log_integration_result(integration_result res) {
	sbmf_log_info("integral: %e", res.integral);
	sbmf_log_info("error: %e", res.error);
	sbmf_log_info("performed evals: %d", res.performed_evals);
	sbmf_log_info("converged: %s", (res.converged) ? "yes" : "no");
}

static void linear_me_integrand(f64* out, f64* in, u32 len, void* data) {
	struct integrand_params* params = data;

	f64 eig1[len];
	f64 eig2[len];
	params->basis.eigenfunc(params->n[0], len, eig1, in);
	params->basis.eigenfunc(params->n[1], len, eig2, in);

	f64 op[len];
	params->op(len, op, in, 0, NULL, params->op_userdata);

	for (u32 i = 0; i < len; ++i) {
		out[i] = eig1[i]*eig2[i]*op[i];
	}
}

static void nonlinear_me_integrand(f64* out, f64* in, u32 len, void* data) {
	struct integrand_params* params = data;

	f64 eig1[len];
	f64 eig2[len];
	params->basis.eigenfunc(params->n[0], len, eig1, in);
	params->basis.eigenfunc(params->n[1], len, eig2, in);

	f64 sample[len*params->component_count];
	for (u32 i = 0; i < params->component_count; ++i) {
		params->basis.sample(params->coeff_count, &params->coeff[i*params->coeff_count], len, &sample[i*len], in);
	}

	f64 op[len];
	params->op(len, op, in, params->component_count, sample, params->op_userdata);

	for (u32 i = 0; i < len; ++i) {
		out[i] = eig1[i]*eig2[i]*op[i];
	}
}

struct gp2c_result gp2c(struct gp2c_settings settings, const u32 component_count, struct gp2c_component component[static component_count]) {
	const u32 N = settings.num_basis_functions;
	const u32 matrix_element_count = N*(N+1)/2;

	/* Setup results struct */
	struct gp2c_result res = {
		.iterations = 0,
		.component_count = component_count,
		.coeff_count = N,
		.coeff = (f64*) sbmf_stack_push(component_count*(N*sizeof(f64))),
		.error = (f64*) sbmf_stack_push(component_count*sizeof(f64)),
		.energy = (f64*) sbmf_stack_push(component_count*sizeof(f64)),
		.hamiltonian = (struct hermitian_bandmat*) sbmf_stack_push(component_count * sizeof(struct hermitian_bandmat)),
	};
	for (u32 i = 0; i < component_count; ++i) {
		res.hamiltonian[i] = hermitian_bandmat_new(N,N);
	}

	/* Place to store coeffs from previous iterations */
	f64* old_coeff = (f64*) sbmf_stack_push(component_count*(N*sizeof(f64)));

	integration_settings int_settings = {
		.gk = (settings.gk.gauss_size > 0) ? settings.gk : gk7,
		.abs_error_tol = 1e-15,
		.rel_error_tol = 1e-15,
		.max_evals = settings.max_iterations,
	};

	/* Setup intial guess values for coeffs */
	for (u32 i = 0; i < component_count; ++i) {
		if (component[i].guess) {
			component[i].guess(&res.coeff[i*res.coeff_count], res.coeff_count);
		} else {
			for (u32 j = 0; j < res.coeff_count; ++j) {
				res.coeff[i*res.coeff_count + j] = 0;
			}
			/* Initialize the i:th component to the i:th eigenfunction */
			res.coeff[i*res.coeff_count + i] = 1;
		}
	}

	struct integrand_params params = {
		.component_count = res.component_count,
		.coeff_count = res.coeff_count,
		.coeff = res.coeff,
		.op = NULL,
		.basis = settings.basis,
	};

	/* Precompute indices for matrix elements
	 * 	This way we get a single array looping
	 * 	over the matrix which is easier to
	 * 	parallelize well.
	 * 		Picking the first hamiltonian as a
	 * 	dummny matrix just to get the correct
	 * 	iteration.
	 * */
	u32 matrix_element_rows[matrix_element_count];
	u32 matrix_element_cols[matrix_element_count];
	{
		u32 matrix_element_index = 0;
		COMPLEX_HERMITIAN_BANDMAT_FOREACH(res.hamiltonian[0], r,c) {
			matrix_element_rows[matrix_element_index] = r;
			matrix_element_cols[matrix_element_index] = c;
			matrix_element_index++;
		}
	}

	/* Construct linear hamiltonian */
	struct hermitian_bandmat linear_hamiltonian = hermitian_bandmat_new_zero(N,N);
	{
		if (settings.ho_potential_perturbation) {
			params.op = settings.ho_potential_perturbation;
#pragma omp parallel for firstprivate(params, int_settings) shared(res)
			for (u32 i = 0; i < matrix_element_count; ++i) {
				u32 r = matrix_element_rows[i];
				u32 c = matrix_element_cols[i];

				int_settings.userdata = &params;
				params.n[0] = r;
				params.n[1] = c;

				integration_result int_res = quadgk_vec(linear_me_integrand, -INFINITY, INFINITY, int_settings);

				//if (fabs(int_res.integral) < 1e-10)
				//	int_res.integral = 0.0;

				if (!int_res.converged) {
					sbmf_log_error("In construction of linear hamiltonian:");
					sbmf_log_error("\tIntegration failed for %d,%d", r,c);
					sbmf_log_integration_result(int_res);
				}
				assert(int_res.converged);

				u32 me_index = hermitian_bandmat_index(linear_hamiltonian, r,c);
				linear_hamiltonian.data[me_index] = int_res.integral;
			}
		}

		for (u32 i = 0; i < res.coeff_count; ++i) {
			u32 me_index = hermitian_bandmat_index(linear_hamiltonian, i,i);
			linear_hamiltonian.data[me_index] += settings.basis.eigenval(i);
		}

		//assert(complex_hermitian_bandmat_is_valid(linear_hamiltonian));
	}

	/* Do the actual iterations */
	for (; res.iterations < settings.max_iterations; ++res.iterations) {
		memcpy(old_coeff, res.coeff, res.component_count * res.coeff_count * sizeof(f64));

		for (u32 i = 0; i < component_count; ++i) {
			params.op = component[i].op;
			params.op_userdata = component[i].userdata;

			/* Construct the i:th component's hamiltonian */
#pragma omp parallel for firstprivate(params, int_settings) shared(res, linear_hamiltonian)
			for (u32 j = 0; j < matrix_element_count; ++j) {
				u32 r = matrix_element_rows[j];
				u32 c = matrix_element_cols[j];
				int_settings.userdata = &params;
				params.n[0] = r;
				params.n[1] = c;

				integration_result int_res = quadgk_vec(nonlinear_me_integrand, -INFINITY, INFINITY, int_settings);

				/* Check if the resultant integral is less than what we can resolve,
				 * if so, zero it.
				 */
				if (fabs(int_res.integral) < 1e-10)
					int_res.integral = 0.0;

				if (!int_res.converged) {
					sbmf_log_error("In construction of component %u's hamiltonian:", i);
					sbmf_log_error("\tIntegration failed for %d,%d", r,c);
					sbmf_log_integration_result(int_res);
				}
				assert(int_res.converged);

				u32 me_index = hermitian_bandmat_index(res.hamiltonian[i], r,c);
				res.hamiltonian[i].data[me_index] = linear_hamiltonian.data[me_index] + int_res.integral;
			}

			assert(hermitian_bandmat_is_valid(res.hamiltonian[i]));
		}

		for (u32 i = 0; i < component_count; ++i) {
			/* Solve for first eigenvector (ground state) */
			struct eigen_result_real eigres = find_eigenpairs_sparse_real(res.hamiltonian[i], 1, EV_SMALLEST_MAG);

			/* Copy energies */
			res.energy[i] = eigres.eigenvalues[0];

			for (u32 j = 0; j < res.coeff_count; ++j) {
				res.coeff[i*res.coeff_count + j] = eigres.eigenvectors[j];
			}
			f64_normalize(&res.coeff[i*res.coeff_count], &res.coeff[i*res.coeff_count], res.coeff_count);
			//f64_normalize(&res.coeff[i*res.coeff_count], &eigres.eigenvectors[0], res.coeff_count);
		}

		/* Calculate error */
		for (u32 i = 0; i < component_count; ++i) {
			f64 sum = 0.0;
			for (u32 j = 0; j < res.coeff_count; ++j) {
				f64 diff = fabs(res.coeff[i*res.coeff_count + j]) - fabs(old_coeff[i*res.coeff_count + j]);
				sum += diff*diff;
			}
			res.error[i] = sqrt(sum);
		}

		/* Break condition */
		sbmf_log_info("gp2c finished iterations %u", res.iterations);
		bool should_exit = true;
		for (u32 i = 0; i < component_count; ++i) {
			if (res.error[i] > settings.error_tol)
				should_exit = false;
			sbmf_log_info("\t[%u] -- error: %.2e, energy: %.2e", i, res.error[i], res.energy[i]);
		}
		if (should_exit)
			break;

#if 1
		/* Call debug callback if requested by user */
		if (settings.measure_every > 0 && res.iterations % settings.measure_every == 0) {
			if (settings.dbgcallback) {
				settings.dbgcallback(settings, res);
			}
		}
#endif
	}

	return res;
}
