#define USE_GSL_INTEGRATION 0

#if USE_GSL_INTEGRATION
	#include <gsl/gsl_integration.h>
#endif

struct integrand_params {
	u32 n[2];
	u32 component_count;
	u32 coeff_count;
	f64* coeff;

	nlse_operator_func* op;
	void* op_userdata;

	struct basis basis;
};

static void sbmf_log_integration_result(struct quadgk_result* res) {
	sbmf_log_info("integral: %.10e", res->integral);
	sbmf_log_info("error: %.10e", res->error);
	//sbmf_log_info("performed evals: %d", res->performed_evals);
	sbmf_log_info("converged: %s", (res->converged) ? "yes" : "no");
}

/* functions to handle spatial guesses */

struct guess_integrand_params {
    u32 n;
	nlse_spatial_guess_func* func;
	struct basis b;
};

void guess_integrand(f64* out, f64* in, u32 len, void* data) {
    struct guess_integrand_params* params = data;

    f64 eig[len];
	params->b.eigenfunc(params->n, len, eig, in);

	f64 sample[len];
	params->func(sample, in, len, NULL);

    for (u32 i = 0; i < len; ++i) {
        out[i] = eig[i] * sample[i];
    }
}

/* function to sample linear and non-linear matrix elements */
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

#if USE_GSL_INTEGRATION
static f64 nonlinear_me_integrand_gsl(f64 in, void* data) {
	struct integrand_params* params = data;

	f64 eig1;
	f64 eig2;
	params->basis.eigenfunc(params->n[0], 1, &eig1, &in);
	params->basis.eigenfunc(params->n[1], 1, &eig2, &in);

	f64 sample[params->component_count];
	for (u32 i = 0; i < params->component_count; ++i) {
		params->basis.sample(params->coeff_count, &params->coeff[i*params->coeff_count], 1, &sample[i], &in);
	}

	f64 op;
	params->op(1, &op, &in, params->component_count, sample, params->op_userdata);

	return eig1*eig2*op;
}
#endif

struct nlse_result nlse_solver(struct nlse_settings settings, const u32 component_count, struct nlse_component* component) {
	/* Lazy */
	const u32 N = settings.num_basis_funcs;

	/* Setup results struct */
	struct nlse_result res = {
		.iterations = 0,
		.component_count = component_count,
		.coeff_count = N,
		.coeff 		 = sbmf_stack_push(component_count*(N*sizeof(f64))),
		.abs_error 	 = sbmf_stack_push(component_count*sizeof(f64)),
		.energy 	 = sbmf_stack_push(component_count*sizeof(f64)),
		.residual 	 = sbmf_stack_push(component_count*sizeof(f64)),
		.hamiltonian = sbmf_stack_push(component_count * sizeof(struct symmetric_bandmat)),
		.converged = false,
	};
	for (u32 i = 0; i < component_count; ++i) {
		res.hamiltonian[i] = symmetric_bandmat_new(N,N);
	}

	/* Place to store coeffs from previous iterations */
	f64* old_coeff = sbmf_stack_push(component_count*(N*sizeof(f64)));
	f64* old_energy = sbmf_stack_push(component_count * sizeof(f64));

	struct quadgk_settings int_settings = {
		.gk = (settings.gk.gauss_size > 0) ? settings.gk : gk15,
		.abs_error_tol = 1e-15,
		.rel_error_tol = 1e-8,
		.max_iters = settings.max_quadgk_iters,
	};

	/* Setup intial guess values for coeffs */
	for (u32 i = 0; i < component_count; ++i) {
		switch (component[i].guess.type) {
			case DEFAULT_GUESS: {
				for (u32 j = 0; j < res.coeff_count; ++j) {
					res.coeff[i*res.coeff_count + j] = 0;
				}
				/* Initialize the i:th component to the i:th eigenfunction */
				res.coeff[i*res.coeff_count + i] = 1;
				break;
			}
			case SPATIAL_GUESS: {
				struct guess_integrand_params p = {
					.func = component[i].guess.data.spatial_guess,
					.b = settings.basis,
				};
				int_settings.userdata = &p;
				f64* out = &res.coeff[i*res.coeff_count];
				for (u32 j = 0; j < res.coeff_count; ++j) {
					p.n = j;
					struct quadgk_result r;

					u8 quadgk_memory[quadgk_required_memory_size(&int_settings)];
					quadgk_infinite_interval(guess_integrand, &int_settings, quadgk_memory, &r);
					out[j] = r.integral;
				}
				f64_normalize(out, out, res.coeff_count);

				break;
			}
			case COEFF_GUESS: {
				component[i].guess.data.coeff_guess(&res.coeff[i*res.coeff_count], res.coeff_count, i);
				break;
			}
			case RANDOM_GUESS: {
				struct symmetric_bandmat bm = symmetric_bandmat_new_zero(N,N);
				SYMMETRIC_BANDMAT_FOREACH(bm, r,c) {
					u32 index = symmetric_bandmat_index(bm, r,c);
					bm.data[index] = 2.0 * ((f64)rand()/(f64)RAND_MAX) - 1.0;
				}
				struct eigen_result_real eigres = find_eigenpairs_sparse_real(bm, 1, EV_SMALLEST);
				for (u32 j = 0; j < res.coeff_count; ++j) {
					res.coeff[i*res.coeff_count + j] = eigres.eigenvectors[j];
				}
				break;
			}
			default:
				sbmf_log_error("nlse_solver: component %u has invalid guess!", i);
				return res;
		};
	}

	struct integrand_params params = {
		.component_count = res.component_count,
		.coeff_count = res.coeff_count,
		.coeff = res.coeff,
		.op = NULL,
		.basis = settings.basis,
	};

	/*
	 * Unique number of me's that need to be calculated in
	 * a symmetric matrix
	 */
	const u32 matrix_element_count = symmetric_bandmat_element_count(N);

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
		SYMMETRIC_BANDMAT_FOREACH(res.hamiltonian[0], r,c) {
			matrix_element_rows[matrix_element_index] = r;
			matrix_element_cols[matrix_element_index] = c;
			matrix_element_index++;
		}
	}

	/* Construct linear hamiltonian */
	struct symmetric_bandmat linear_hamiltonian = symmetric_bandmat_new_zero(N,N);
	{
		if (settings.spatial_pot_perturbation) {
			params.op = settings.spatial_pot_perturbation;
#pragma omp parallel for firstprivate(params, int_settings) shared(res)
			for (u32 i = 0; i < matrix_element_count; ++i) {
				u32 r = matrix_element_rows[i];
				u32 c = matrix_element_cols[i];

				int_settings.userdata = &params;
				params.n[0] = r;
				params.n[1] = c;

				u8 quadgk_memory[quadgk_required_memory_size(&int_settings)];
				struct quadgk_result int_res;
				quadgk_infinite_interval(linear_me_integrand, &int_settings, quadgk_memory, &int_res);

				if (fabs(int_res.integral) <= settings.zero_threshold)
					int_res.integral = 0.0;

				if (!int_res.converged) {
					sbmf_log_error("In construction of linear hamiltonian:");
					sbmf_log_error("\tIntegration failed for %d,%d", r,c);
					sbmf_log_integration_result(&int_res);
				}
				assert(int_res.converged);

				u32 me_index = symmetric_bandmat_index(linear_hamiltonian, r,c);
				linear_hamiltonian.data[me_index] = int_res.integral;
			}
		}

		for (u32 i = 0; i < res.coeff_count; ++i) {
			u32 me_index = symmetric_bandmat_index(linear_hamiltonian, i,i);
			linear_hamiltonian.data[me_index] += settings.basis.eigenval(i);
		}

		assert(symmetric_bandmat_is_valid(linear_hamiltonian));
	}

#if USE_GSL_INTEGRATION
	const u32 num_threads = omp_get_max_threads();
	gsl_integration_workspace* ws[num_threads];
	/* Create gsl workspaces if necessary */
	for (u32 i = 0; i < num_threads; ++i) {
		ws[i] = gsl_integration_workspace_alloc(settings.max_iterations);
	}
#endif

	struct symmetric_bandmat old_Hs[res.component_count];
	for (u32 i = 0; i < res.component_count; ++i) {
		old_Hs[i] = symmetric_bandmat_new(N,N);
	}

	/* Do the actual iterations */
	for (; res.iterations < settings.max_iterations; ++res.iterations) {
		sbmf_log_info("nlse starting iteration %u", res.iterations);
		memcpy(old_coeff, res.coeff, res.component_count * N * sizeof(f64));
		memcpy(old_energy, res.energy, res.component_count * sizeof(f64));

		/* Call debug callback if requested by user */
		if (settings.measure_every > 0 &&
				settings.debug_callback &&
				res.iterations % settings.measure_every == 0) {
			settings.debug_callback(settings, res);
		}

		/* Should we still apply mixing? */
		if (settings.orbital_mixing > 0 && settings.mix_until_iteration > 0 && res.iterations == settings.mix_until_iteration) {
			settings.orbital_mixing = 0;
		}
		if (settings.hamiltonian_mixing > 0 && settings.mix_until_iteration > 0 && res.iterations == settings.mix_until_iteration) {
			settings.hamiltonian_mixing = 0;
		}

		if (settings.orbital_choice != NLSE_ORBITAL_MAXIMUM_OVERLAP && settings.mom_enable_at_iteration > 0 && res.iterations == settings.mom_enable_at_iteration) {
			settings.orbital_choice = NLSE_ORBITAL_MAXIMUM_OVERLAP;
		}

		/*
		 * Construct all the Hamiltonians
		 */
		for (u32 i = 0; i < component_count; ++i) {
			params.op = component[i].op;
			params.op_userdata = component[i].userdata;

			memcpy(old_Hs[i].data, res.hamiltonian[i].data, N*N*sizeof(f64));

			/* Construct the i:th component's hamiltonian */
#pragma omp parallel for firstprivate(params, int_settings) shared(res, linear_hamiltonian)
			for (u32 j = 0; j < matrix_element_count; ++j) {
				u32 r = matrix_element_rows[j];
				u32 c = matrix_element_cols[j];
				int_settings.userdata = &params;
				params.n[0] = r;
				params.n[1] = c;

#if USE_GSL_INTEGRATION
				integration_result int_res;
				gsl_function F;
				F.function = nonlinear_me_integrand_gsl;
				F.params = &params;
				gsl_integration_qagi(&F, 1e-10, 1e-7, settings.max_iterations, ws[omp_get_thread_num()], &int_res.integral, &int_res.error);
				int_res.converged = true;
#else
				u8 quadgk_memory[quadgk_required_memory_size(&int_settings)];

				struct quadgk_result int_res;
				quadgk_infinite_interval(nonlinear_me_integrand, &int_settings, quadgk_memory, &int_res);
#endif

				/* Check if the resultant integral is less than what we can resolve,
				 * if so, zero it.
				 */
				if (fabs(int_res.integral) <= settings.zero_threshold)
					int_res.integral = 0.0;

				if (!int_res.converged) {
					sbmf_log_error("In construction of component %u's hamiltonian:", i);
					sbmf_log_error("\tIntegration failed for %d,%d", r,c);
					sbmf_log_integration_result(&int_res);
				}
				assert(int_res.converged);

				u32 me_index = symmetric_bandmat_index(res.hamiltonian[i], r,c);
				res.hamiltonian[i].data[me_index] = linear_hamiltonian.data[me_index] + int_res.integral;

				if (settings.hamiltonian_mixing > 0) {
					res.hamiltonian[i].data[me_index] = (1.0 - settings.hamiltonian_mixing)*res.hamiltonian[i].data[me_index] + settings.hamiltonian_mixing*old_Hs[i].data[me_index];
				}
			}

			assert(symmetric_bandmat_is_valid(res.hamiltonian[i]));
		}

		/*
		 * Solve and normalize all the Hamiltonian
		 * eigenvalue problems
		 */
		for (u32 i = 0; i < component_count; ++i) {
			struct eigen_result_real eigres;
			u32 new_orbital_index = 0;
			switch (settings.orbital_choice) {
				case NLSE_ORBITAL_LOWEST_ENERGY: {
					eigres = find_eigenpairs_sparse_real(res.hamiltonian[i], 1, EV_SMALLEST);
					f64_normalize(&eigres.eigenvectors[0], &eigres.eigenvectors[0], res.coeff_count);
					break;
				}
				case NLSE_ORBITAL_MAXIMUM_OVERLAP: {
					eigres = find_eigenpairs_sparse_real(res.hamiltonian[i], settings.mom_orbitals_to_consider, EV_SMALLEST);
					f64 maximum_overlap = -INFINITY;
					for (u32 j = 0; j < settings.mom_orbitals_to_consider; ++j) {
						f64_normalize(&eigres.eigenvectors[j*res.coeff_count], &eigres.eigenvectors[j*res.coeff_count], res.coeff_count);
						f64 overlap = 0.0;
						for (u32 k = 0; k < res.coeff_count; ++k) {
							f64 a = eigres.eigenvectors[j*res.coeff_count + k];
							f64 b = old_coeff[i*res.coeff_count + k];
							overlap += a*b;
						}
						overlap = fabs(overlap);

						if (overlap > maximum_overlap) {
							new_orbital_index = j;
							maximum_overlap = overlap;
						}
					}
					break;
				}
				default: {
					sbmf_log_panic("nlse_solver(...): Unspecified energy selection method!");
					assert(0);
				}
			};

			/* Copy energies and coeffs */
			res.energy[i] = eigres.eigenvalues[new_orbital_index];
			for (u32 j = 0; j < res.coeff_count; ++j) {
				res.coeff[i*res.coeff_count + j] = eigres.eigenvectors[res.coeff_count*new_orbital_index + j];
			}
			f64_normalize(&res.coeff[i*res.coeff_count], &res.coeff[i*res.coeff_count], res.coeff_count);

			/* Apply orbital mixing */
			if (settings.orbital_mixing > 0.0) {
				res.energy[i] = (1.0 - settings.orbital_mixing) * res.energy[i] + settings.orbital_mixing*old_energy[i];

				for (u32 j = 0; j < res.coeff_count; ++j) {
					res.coeff[i*res.coeff_count + j] =
						(1.0 - settings.orbital_mixing) * res.coeff[i*res.coeff_count + j] + settings.orbital_mixing * old_coeff[i*res.coeff_count + j];
				}
				f64_normalize(&res.coeff[i*res.coeff_count], &res.coeff[i*res.coeff_count], res.coeff_count);
			}

		}

		/* Calculate error */
		for (u32 i = 0; i < component_count; ++i) {
			f64 sum_diff = 0.0;
			f64 sum_prev = 0.0;
			for (u32 j = 0; j < res.coeff_count; ++j) {
				f64 diff = fabs(res.coeff[i*res.coeff_count + j]) - fabs(old_coeff[i*res.coeff_count + j]);
				sum_diff += diff*diff;
				sum_prev += old_coeff[i*res.coeff_count + j]*old_coeff[i*res.coeff_count + j];
			}
			res.abs_error[i] = sqrt(sum_diff);
		}

		/* Break condition */
		bool should_exit = true;
		for (u32 i = 0; i < component_count; ++i) {
			if (res.abs_error[i] > settings.abs_error_tol)
				should_exit = false;
			sbmf_log_info("\t[%u] -- abs error: %.10e energy: %.10e", i, res.abs_error[i], res.energy[i]);
		}

		if (should_exit) {
			res.converged = true;
			break;
		}
	}

#if USE_GSL_INTEGRATION
	/* free gsl workspaces if necessary */
	for (u32 i = 0; i < num_threads; ++i) {
		gsl_integration_workspace_free(ws[i]);
	}
#endif

	/* Compute residuals */
	for (u32 i = 0; i < component_count; ++i) {
		f64 ans1[N], ans2[N];
		symmetric_bandmat_mulv(ans1, res.hamiltonian[i], &res.coeff[i*N]);

		for (u32 j = 0; j < N; ++j) {
			ans2[j] = res.energy[i] * res.coeff[i*N + j];
		}

		f64 sum = 0.0;
		for (u32 j = 0; j < N; ++j) {
			sum += fabs(ans1[j] - ans2[j]);
		}
		sum = sqrt(sum);

		res.residual[i] = sum;
		sbmf_log_info("\t[%u] residual %e", i, sum);
	}

	return res;
}

/*
 * basic serialization
 */

void nlse_write_to_binary_file(const char* file, struct nlse_result res) {
	FILE* fd = fopen(file, "w");
	fwrite(&res.iterations, 	 sizeof(u32), 1, fd);
	fwrite(&res.component_count, sizeof(u32), 1, fd);
	fwrite(&res.coeff_count, 	 sizeof(u32), 1, fd);

	fwrite(res.coeff,  sizeof(f64), res.coeff_count*res.component_count, fd);
	fwrite(res.abs_error,  sizeof(f64), res.component_count, fd);
	fwrite(res.energy, sizeof(f64), res.component_count, fd);

	for (u32 i = 0; i < res.component_count; ++i) {
		fwrite(&res.hamiltonian[i].size, 	  sizeof(u32), 1, fd);
		fwrite(&res.hamiltonian[i].bandcount, sizeof(u32), 1, fd);
		u32 elements_written = fwrite(res.hamiltonian[i].data, sizeof(f64),
				res.hamiltonian[i].bandcount*res.hamiltonian[i].size,
				fd);
		assert(elements_written == res.hamiltonian[i].bandcount*res.hamiltonian[i].size);
	}
	fwrite(&res.converged, sizeof(bool), 1, fd);
	fclose(fd);
}

struct nlse_result nlse_read_from_binary_file(const char* file) {
	struct nlse_result res;

	FILE* fd = fopen(file, "r");
	fread(&res.iterations, 	 	sizeof(u32), 1, fd);
	fread(&res.component_count, sizeof(u32), 1, fd);
	fread(&res.coeff_count,  	sizeof(u32), 1, fd);

	res.coeff = sbmf_stack_push(sizeof(f64)*res.coeff_count*res.component_count);
	fread(res.coeff,  sizeof(f64), res.coeff_count*res.component_count, fd);

	res.abs_error = sbmf_stack_push(sizeof(f64)*res.component_count);
	fread(res.abs_error,  sizeof(f64), res.component_count, fd);

	res.energy = sbmf_stack_push(sizeof(f64)*res.component_count);
	fread(res.energy, sizeof(f64), res.component_count, fd);

	res.hamiltonian = sbmf_stack_push(res.component_count * sizeof(struct symmetric_bandmat));
	for (u32 i = 0; i < res.component_count; ++i) {
		u32 size, bandcount;
		fread(&size, 	  sizeof(u32), 1, fd);
		fread(&bandcount, sizeof(u32), 1, fd);
		res.hamiltonian[i] = symmetric_bandmat_new(bandcount, size);
		u32 bytes_written = fread(res.hamiltonian[i].data, sizeof(f64),
				bandcount*size,
				fd);
		assert(bytes_written == bandcount*size);
	}
	fread(&res.converged, sizeof(bool), 1, fd);
	fclose(fd);

	sbmf_log_info("loaded:");
	sbmf_log_info("    iterations: %u", res.iterations);
	sbmf_log_info("    component_count: %u", res.component_count);
	sbmf_log_info("    coeff_count: %u", res.coeff_count);
	for (u32 i = 0; i < res.component_count; ++i) {
		sbmf_log_info("    [%u]: size: %u", i, res.hamiltonian[i].size);
		sbmf_log_info("    [%u]: bandcount: %u", i, res.hamiltonian[i].bandcount);
	}

	return res;
}
