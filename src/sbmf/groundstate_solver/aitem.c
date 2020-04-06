#include "groundstate_solver.h"
#include <sbmf/common/math_utils.h>

#include <fftw3.h>

#include <string.h>
#include <stdlib.h>

gss_result aitem_execute(gss_settings settings, gss_potential_func* potential, gss_guess_func* guess) {
	const real_t dx = settings.length_x / settings.resolution;
	const real_t dy = settings.length_y / settings.resolution;
	const int_t N = settings.resolution;
	// applying fft -> ifft in fftw scales input by N = n0*n1*...*nk.
	// this is needed to cancel that scaling.
	const real_t ifft_factor = 1.0/(N*N);

	gss_result result = {
		.settings = settings,
		.rows = settings.resolution,
		.cols = settings.resolution,
		.wavefunction = malloc(N*N*sizeof(complex_t)),
		.error = 0.0,
		.iterations = 0,
	};

	real_t X[N*N];
	real_t Y[N*N];
	real_t KX[N*N];
	real_t KY[N*N];
	{
		FOREACH_ROW(N,N, i,j) {
			X[mat_idx(N,i,j)] = -settings.length_x/2 + j*dx;
			Y[mat_idx(N,i,j)] = -settings.length_y/2 + i*dy;
			KX[mat_idx(N,i,j)] = (j<N/2) ? 2*M_PI/settings.length_x*j : 2*M_PI/settings.length_x*(-N+j);
			KY[mat_idx(N,i,j)] = (i<N/2) ? 2*M_PI/settings.length_y*i : 2*M_PI/settings.length_y*(-N+i);
		}
	}

	// Variables used in computation
	//fftw_complex U[N*N];
	fftw_complex old_wavefunction[N*N];
	fftw_complex MinvU[N*N];
	fftw_complex L00U[N*N];
	fftw_complex fftU[N*N];

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
			real_t x = X[mat_idx(N,i,j)];
			real_t y = Y[mat_idx(N,i,j)];
			result.wavefunction[mat_idx(N,i,j)] = guess(x,y);
		}

		// Normalize w.f.
		real_t d = settings.A/max_abs_complex(N,N,result.wavefunction);
		FOREACH_ROW(N,N, i,j) {
			result.wavefunction[mat_idx(N,i,j)] *= d;
		}
	}

	for (; result.iterations < settings.max_iterations; ++result.iterations)  {
		memcpy(old_wavefunction, result.wavefunction, N*N*sizeof(fftw_complex));
		
		// Compute L00U
		{
			// Compute elementwise term
			FOREACH_ROW(N,N, i,j) {
				L00U[mat_idx(N,i,j)] = potential(X[mat_idx(N,i,j)], Y[mat_idx(N,i,j)], result.wavefunction[mat_idx(N,i,j)]);
			}

			// Compute fftU
			memcpy(fft_in, result.wavefunction, N*N*sizeof(complex_t));
			fftw_execute(fft_plan);
			memcpy(fftU, fft_out, N*N*sizeof(complex_t));

			// Compute and add the fft term
			FOREACH_ROW(N,N, i,j) {
				real_t kx = KX[mat_idx(N,i,j)];
				real_t ky = KY[mat_idx(N,i,j)];
				ifft_in[mat_idx(N,i,j)] = -(kx*kx + ky*ky)*fftU[mat_idx(N,i,j)];
			}
			fftw_execute(ifft_plan);
			FOREACH_ROW(N,N, i,j) {
				L00U[mat_idx(N,i,j)] += ifft_factor*ifft_out[mat_idx(N,i,j)];
			}
		}

		// Compute MinvU
		{
			FOREACH_ROW(N,N, i,j) {
				real_t kx = KX[mat_idx(N,i,j)];
				real_t ky = KY[mat_idx(N,i,j)];
				real_t factor = settings.c + kx*kx + ky*ky;
				ifft_in[mat_idx(N,i,j)] = fftU[mat_idx(N,i,j)] / factor;
			}
			fftw_execute(ifft_plan);
			FOREACH_ROW(N,N, i,j) {
				MinvU[mat_idx(N,i,j)] = ifft_factor*ifft_out[mat_idx(N,i,j)];
			}
		}
		
		// Update and normalize U
		{
			fftw_complex mu_denom = 0.0;
			fftw_complex mu_numer = 0.0;
			FOREACH_ROW(N,N, i,j) {
				mu_numer += (L00U[mat_idx(N,i,j)]*MinvU[mat_idx(N,i,j)]);
				mu_denom += (result.wavefunction[mat_idx(N,i,j)]*MinvU[mat_idx(N,i,j)]);
			}
			fftw_complex mu = mu_numer/mu_denom;

			FOREACH_ROW(N,N, i,j) {
				fft_in[mat_idx(N,i,j)] = L00U[mat_idx(N,i,j)] - mu*result.wavefunction[mat_idx(N,i,j)];
			}
			fftw_execute(fft_plan);
			FOREACH_ROW(N,N, i,j) {
				real_t kx = KX[mat_idx(N,i,j)];
				real_t ky = KY[mat_idx(N,i,j)];
				real_t factor = settings.c + kx*kx + ky*ky;
				ifft_in[mat_idx(N,i,j)] = fft_out[mat_idx(N,i,j)] / factor;
			}
			fftw_execute(ifft_plan);

			FOREACH_ROW(N,N, i,j) {
				result.wavefunction[mat_idx(N,i,j)] += ifft_factor*ifft_out[mat_idx(N,i,j)]*settings.dt;
			}

			// Normalize
			real_t d = settings.A/max_abs_complex(N,N,result.wavefunction);
			FOREACH_ROW(N,N, i,j) {
				result.wavefunction[mat_idx(N,i,j)] *= d;
			}
		}

		// Calculate error
		{
			double sum = 0.0;
			FOREACH_ROW(N,N, i,j) {
				double diff = cabs(result.wavefunction[mat_idx(N,i,j)]-old_wavefunction[mat_idx(N,i,j)]);
				sum += diff*diff;
			}
			result.error = sqrt(sum*dx*dy);

			if (result.error < settings.error_tol) {
				break;
			}
		}
	}

	// Cleanup
	fftw_destroy_plan(fft_plan);
	fftw_destroy_plan(ifft_plan);
	fftw_free(fft_in);
	fftw_free(fft_out);
	fftw_free(ifft_in);
	fftw_free(ifft_out);

	return result;
}
