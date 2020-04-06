#include "groundstate_solver.h"
#include <sbmf/common/math_utils.h>

#include <fftw3.h>

#include <string.h>
#include <stdlib.h>

static real_t VK(real_t kx, real_t ky, complex_t u) {
	return 0.5*(kx*kx + ky*ky);
}

static inline void apply_step_op(real_t ds, real_t dt, complex_t* out, gss_potential_func* potential, int_t rows, int_t cols, real_t* xmat, real_t* ymat, complex_t* wfmat) {
	FOREACH_ROW(rows,cols, i,j) {
		int_t idx = mat_idx(rows,i,j);
		real_t potval = potential(xmat[idx], ymat[idx], wfmat[idx]);
		real_t opval  = exp(-potval*dt);
		out[idx] = wfmat[idx] * opval * ds;
	}
}

gss_result item_execute(gss_settings settings, gss_potential_func* potential, gss_guess_func* guess) {
	const real_t Lx = settings.length_x;
	const real_t Ly = settings.length_y;
	const real_t dx = settings.length_x / settings.resolution;
	const real_t dy = settings.length_y / settings.resolution;
	const real_t dt = settings.dt;
	const real_t ds = (dx*dy);
	const int_t N = settings.resolution;
	// applying fft -> ifft in fftw scales input by N = n0*n1*...*nk.
	// this is needed to cancel that scaling.
	const real_t ifft_factor = 1.0/(N*N);

	gss_debug dbg;
	if (settings.measure_every > 0) {
		dbg.count = settings.max_iterations/settings.measure_every;
		dbg.current = 0;
		dbg.error = malloc(dbg.count * sizeof(real_t));
		dbg.iteration_time = malloc(dbg.count * sizeof(real_t));
		dbg.wavefunction = malloc((N*N) * dbg.count * sizeof(complex_t));
	}

	gss_result result = {
		.settings = settings,
		.debug = dbg,
		.rows = settings.resolution,
		.cols = settings.resolution,
		.X = malloc(N*N*sizeof(real_t)),
		.Y = malloc(N*N*sizeof(real_t)),
		.wavefunction = malloc(N*N*sizeof(complex_t)),
		.error = 0.0,
		.iterations = 0,
	};

	real_t* X = result.X;
	real_t* Y = result.Y;
	real_t KX[N*N];
	real_t KY[N*N];
	FOREACH_ROW(N,N, i,j) {
		int_t idx = mat_idx(N,i,j);
		X[idx] = -settings.length_x/2 + j*dx;
		Y[idx] = -settings.length_y/2 + i*dy;
		KX[idx] = (j<N/2) ? 2*M_PI/settings.length_x*j : 2*M_PI/settings.length_x*(-N+j);
		KY[idx] = (i<N/2) ? 2*M_PI/settings.length_y*i : 2*M_PI/settings.length_y*(-N+i);
	}

	// Variables used in computation
	//fftw_complex U[N*N];
	fftw_complex old_wavefunction[N*N];

	// FFT def.
	complex_t* fft_in   = fftw_malloc(sizeof(complex_t) * N*N);
	complex_t* fft_out  = fftw_malloc(sizeof(complex_t) * N*N);
	complex_t* ifft_in  = fftw_malloc(sizeof(complex_t) * N*N);
	complex_t* ifft_out = fftw_malloc(sizeof(complex_t) * N*N);

	fftw_plan fft_plan  = fftw_plan_dft_2d(N, N, fft_in, fft_out, FFTW_FORWARD, FFTW_MEASURE);
	fftw_plan ifft_plan = fftw_plan_dft_2d(N, N, ifft_in, ifft_out, FFTW_BACKWARD, FFTW_MEASURE);

	// W.f. def.
	{
		FOREACH_ROW(N,N, i,j) {
			int_t idx = mat_idx(N,i,j);
			real_t x = X[idx];
			real_t y = Y[idx];
			result.wavefunction[idx] = guess(x,y);
		}

		real_t sum = 0.0;
		FOREACH_ROW(N,N, i,j) {
			int_t idx = mat_idx(N,i,j);
			real_t tmp = cabs(result.wavefunction[idx]);
			sum += tmp*tmp;
		}
		real_t scaling = 1.0/sqrt(sum*ds);
		FOREACH_ROW(N,N, i,j) {
			int_t idx = mat_idx(N,i,j);
			result.wavefunction[idx] *= scaling;
		}
	}

	apply_step_op(1.0, dt/2.0, fft_in, potential, N,N, X,Y, result.wavefunction);
	fftw_execute(fft_plan);

	apply_step_op(ds, dt, ifft_in, VK, N,N, KX, KY, fft_out);
	fftw_execute(ifft_plan);

	FOREACH_ROW(N,N, i,j) {
		int_t idx = mat_idx(N,i,j);
		result.wavefunction[idx] = ifft_factor*ifft_out[idx]/(Lx*Ly);
	}

	for (; result.iterations < settings.max_iterations; ++result.iterations)  {
		memcpy(old_wavefunction, result.wavefunction, N*N*sizeof(fftw_complex));

		// Apply operators
		{
			apply_step_op(1.0, dt, fft_in, potential, N,N, X,Y, result.wavefunction);
			fftw_execute(fft_plan);

			apply_step_op(ds, dt, ifft_in, VK, N,N, KX, KY, fft_out);
			fftw_execute(ifft_plan);

			FOREACH_ROW(N,N, i,j) {
				int_t idx = mat_idx(N,i,j);
				result.wavefunction[idx] = ifft_factor*ifft_out[idx]/(Lx*Ly);
			}
		}
		
		// Normalize wavefunction
		{
			real_t sum = 0.0;
			FOREACH_ROW(N,N, i,j) {
				int_t idx = mat_idx(N,i,j);
				real_t tmp = cabs(result.wavefunction[idx]);
				sum += tmp*tmp;
			}
			real_t scaling = 1.0/sqrt(sum*ds);
			FOREACH_ROW(N,N, i,j) {
				int_t idx = mat_idx(N,i,j);
				result.wavefunction[idx] *= scaling;
			}
		}

		// Calculate error
		{
			real_t sum = 0.0;
			FOREACH_ROW(N,N, i,j) {
				int_t idx = mat_idx(N,i,j);
				real_t diff = cabs(result.wavefunction[idx]-old_wavefunction[idx]);
				sum += diff*diff;
			}
			result.error = sqrt(sum*ds);
		}

		if (settings.measure_every > 0 && result.iterations % settings.measure_every == 0) {
			result.debug.error[result.debug.current] = result.error;
			memcpy(&result.debug.wavefunction[result.debug.current*(N*N)], result.wavefunction, N*N*sizeof(complex_t));
			result.debug.current++;
		}

		if (result.error < settings.error_tol) {
			break;
		}
	}

	apply_step_op(1.0, dt/2.0, fft_in, potential, N,N, X,Y, result.wavefunction);
	
	// Cleanup
	fftw_destroy_plan(fft_plan);
	fftw_destroy_plan(ifft_plan);
	fftw_free(fft_in);
	fftw_free(fft_out);
	fftw_free(ifft_in);
	fftw_free(ifft_out);

	return result;
}
