#include "groundstate_solver.h"
#include <sbmf/common/math_utils.h>
#include <sbmf/common/fdm.h>
#include <lapacke.h>
#include <assert.h>
#include <string.h>

gss_result niter_execute(gss_settings settings,
													 gss_potential_func* potential, gss_guess_func* guess) {
	const real_t dx = settings.length_x / settings.resolution;
	const real_t dy = settings.length_y / settings.resolution;
	const int_t N = settings.resolution;

	gss_result result = {
		.settings = settings,
		.rows = settings.resolution,
		.cols = settings.resolution,
		.wavefunction = malloc(N* N*sizeof(complex_t)),
		.error = 0.0,
		.iterations = 0,
	};

	real_t X[N*N];
	real_t Y[N*N];
	FOREACH_ROW(N,N, i,j) {
		X[mat_idx(N,i,j)] = -settings.length_x/2 + j*dx;
		Y[mat_idx(N,i,j)] = -settings.length_y/2 + i*dy;
	}

	FOREACH_ROW(N,N, i,j) {
		real_t x = X[mat_idx(N,i,j)];
		real_t y = Y[mat_idx(N,i,j)];
		result.wavefunction[mat_idx(N,i,j)] = guess(x,y);
	}

	bandmat finitediffmat = generate_fd_matrix(N, -4, 2, (real_t[]){dx, dy});

	complex_t* old_wavefunction = malloc(N*N*sizeof(complex_t));
	complex_t* temp_fdmat = malloc(
														finitediffmat.bandcount*finitediffmat.size*sizeof(complex_t));

	real_t* outD = malloc(finitediffmat.size*sizeof(real_t));
	real_t* outE = malloc((finitediffmat.size-1) * sizeof(real_t));
	complex_t* outQ = malloc(finitediffmat.order*sizeof(complex_t));

	{
		fclose(fopen("aitem_wf", "w"));
		fclose(fopen("aitem_error_output", "w"));
	}
	FILE* fdwf = fopen("aitem_wf", "a");
	FILE* fder = fopen("aitem_error_output", "a");

	for (; result.iterations < settings.max_iterations; ++result.iterations)  {
		if (result.iterations % 10 == 0)
			printf("iteration %d/%d [err: %e]\n", result.iterations, settings.max_iterations, result.error);

		memcpy(old_wavefunction, result.wavefunction, N*N*sizeof(complex_t));
		memcpy(temp_fdmat, finitediffmat.bands,
					 finitediffmat.bandcount*finitediffmat.size*sizeof(complex_t));

		for (int_t i = 0; i < finitediffmat.size; ++i) {
			temp_fdmat[i + (finitediffmat.bandcount-1)*finitediffmat.size] += potential(
						X[i], Y[i], result.wavefunction[i]);
		}

		{
			int res = 0;

			// Reduce to tridiag-form
			res = LAPACKE_zhbtrd(LAPACK_ROW_MAJOR, 'V', 'U', finitediffmat.size,
													 finitediffmat.bandcount-1, temp_fdmat, finitediffmat.size, outD, outE,
													 outQ, finitediffmat.size);
			assert(res == 0);

			// Solve eigenvalue problem
			res = LAPACKE_zsteqr(LAPACK_ROW_MAJOR, 'V', finitediffmat.size, outD, outE,
													 outQ, finitediffmat.size);
			assert(res == 0);
		}

		// Eigenvalue in outD
		// Eigenvectors in outQ

		// Take first wavefunc corresponding to lowest energy and use that one for next iter.
		for (int_t i = 0; i < finitediffmat.size; ++i) {
			result.wavefunction[i] = outQ[0 + i*finitediffmat.size];
		}

		// Compute error
		{
			double sum = 0.0;
			FOREACH_ROW(N,N, i,j) {
				double diff = cabs(result.wavefunction[mat_idx(N,i,j)]-old_wavefunction[mat_idx(N,i,
													 j)]);
				sum += diff*diff;
			}
			result.error = sqrt(sum*dx*dy);

			if (result.error < settings.error_tol) {
				break;
			}
		}

		{
			for (int_t i = 0; i < N*N; ++i) {
				real_t cr = creal(result.wavefunction[i]);
				real_t ci = cimag(result.wavefunction[i]);
				fprintf(fdwf, "%lf%s%lfi",
						cr,
						(ci < 0.0) ? "-" : "+",
						fabs(ci));

				if (i < N*N-1)
					fprintf(fdwf, "\t");
				else
					fprintf(fdwf, "\n");
			}
		}
		{
			fprintf(fder, "%lf\n", result.error);
		}
	}

	fclose(fdwf);
	fclose(fder);

	free(temp_fdmat);
	free(old_wavefunction);
	free(finitediffmat.bands);

	return result;
}
