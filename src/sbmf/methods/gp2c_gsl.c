#include <sbmf/methods/gp2c.h>
#include <sbmf/sbmf.h>
#include <sbmf/math/find_eigenpairs.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/math/functions.h>

#include <gsl/gsl_integration.h>
#include <omp.h>

#include <assert.h> /* Not the correct way to handle this */

struct integrand_params {
	u32 n[2];
	u32 component_count;
	u32 coeff_count;
	c64* coeff;
	gp2c_operator_func* op;
	struct basis basis;
};


static f64 linear_me_integrand(f64 x, void* data) {
	struct integrand_params* params = data;

	f64 eig1;
	f64 eig2;
	params->basis.eigenfunc(params->n[0], 1, &eig1, &x);
	params->basis.eigenfunc(params->n[1], 1, &eig2, &x);

	f64 op;
	params->op(1, &op, &x, 0, NULL);

	return eig1*eig2*op;
}

static f64 nonlinear_me_integrand(f64 x, void* data) {
	struct integrand_params* params = data;

	f64 eig1;
	f64 eig2;
	params->basis.eigenfunc(params->n[0], 1, &eig1, &x);
	params->basis.eigenfunc(params->n[1], 1, &eig2, &x);

	c64 sample[params->component_count];
	for (u32 i = 0; i < params->component_count; ++i) {
		params->basis.sample(params->coeff_count, &params->coeff[i*params->coeff_count], 1, &sample[i], &x);
	}

	f64 op;
	params->op(1, &op, &x, params->component_count, sample);

	return eig1*eig2*op;
}






struct gp2c_result gp2c_gsl(struct gp2c_settings settings, const u32 component_count, struct gp2c_component component[static component_count]) {
	const u32 N = settings.num_basis_functions;
	const u32 matrix_element_count = N*(N+1)/2;

	/* Setup results struct */
	struct gp2c_result res = {
		.iterations = 0,
		.component_count = component_count,
		.coeff_count = N,
		.coeff = (c64*) sbmf_stack_push(component_count*(N*sizeof(c64))),
		.error = (f64*) sbmf_stack_push(component_count*sizeof(f64)),
		.energy = (f64*) sbmf_stack_push(component_count*sizeof(f64)),
		.hamiltonian = (struct hermitian_bandmat*) sbmf_stack_push(component_count * sizeof(struct hermitian_bandmat)),
	};
	for (u32 i = 0; i < component_count; ++i) {
		res.hamiltonian[i] = hermitian_bandmat_new(N,N);
	}

	/* Place to store coeffs from previous iterations */
	c64* old_coeff = (c64*) sbmf_stack_push(component_count*(N*sizeof(c64)));

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

	u32 limit = 1e8;
	f64 atol = 1e-15;
	f64 rtol = 1e-7;

	/* Construct thread data */
	gsl_integration_workspace* ws[omp_get_max_threads()];
	for (u32 i = 0; i < omp_get_max_threads(); ++i) {
		ws[i] = gsl_integration_workspace_alloc(limit);
	}

	/* Construct linear hamiltonian */
	struct hermitian_bandmat linear_hamiltonian = hermitian_bandmat_new_zero(N,N);
	{
		if (settings.ho_potential_perturbation) {
			params.op = settings.ho_potential_perturbation;
#pragma omp parallel for firstprivate(params) shared(res)
			for (u32 i = 0; i < matrix_element_count; ++i) {
				u32 r = matrix_element_rows[i];
				u32 c = matrix_element_cols[i];

				params.n[0] = r;
				params.n[1] = c;

				gsl_function F;
				F.function = linear_me_integrand;
				F.params = &params;

				f64 result;
				f64 error;
				gsl_integration_qagi(&F, atol, rtol, limit, ws[omp_get_thread_num()], &result, &error);

				u32 me_index = hermitian_bandmat_index(linear_hamiltonian, r,c);
				linear_hamiltonian.data[me_index] = result;
			}
		}

		for (u32 i = 0; i < res.coeff_count; ++i) {
			u32 me_index = hermitian_bandmat_index(linear_hamiltonian, i,i);
			linear_hamiltonian.data[me_index] += settings.basis.eigenval(i);
		}
	}

	/* Do the actual iterations */
	for (; res.iterations < settings.max_iterations; ++res.iterations) {
		memcpy(old_coeff, res.coeff, res.component_count * res.coeff_count * sizeof(c64));

		for (u32 i = 0; i < component_count; ++i) {
			params.op = component[i].op;

			/* Construct the i:th component's hamiltonian */
#pragma omp parallel for firstprivate(params) shared(res, linear_hamiltonian)
			for (u32 j = 0; j < matrix_element_count; ++j) {
				u32 r = matrix_element_rows[j];
				u32 c = matrix_element_cols[j];

				params.n[0] = r;
				params.n[1] = c;

				gsl_function F;
				F.function = nonlinear_me_integrand;
				F.params = &params;

				f64 result;
				f64 error;
				gsl_integration_qagi(&F, atol, rtol, limit, ws[omp_get_thread_num()], &result, &error);

				if (fabs(result) <= error || fabs(result) <= 1e-10) {
					result = 0.0;
				}

				u32 me_index = hermitian_bandmat_index(res.hamiltonian[i], r,c);
				res.hamiltonian[i].data[me_index] = linear_hamiltonian.data[me_index] + result;
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
			c64_normalize(&res.coeff[i*res.coeff_count], &res.coeff[i*res.coeff_count], res.coeff_count);
			//f64_normalize(&res.coeff[i*res.coeff_count], &eigres.eigenvectors[0], res.coeff_count);
		}

		/* Calculate error */
		for (u32 i = 0; i < component_count; ++i) {
			f64 sum = 0.0;
			for (u32 j = 0; j < res.coeff_count; ++j) {
				f64 diff = cabs(res.coeff[i*res.coeff_count + j]) - cabs(old_coeff[i*res.coeff_count + j]);
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

	/* Free gsl stuff */
	for (u32 i = 0; i < omp_get_max_threads(); ++i) {
		gsl_integration_workspace_free(ws[i]);
	}

	return res;
}
