#include "groundstate_solver.h"
#include <sbmf/common/math_utils.h>

#include <fftw3.h>

#include <string.h>
#include <stdlib.h>
#include <assert.h>

static real_t VK(real_t* v, int_t n, complex_t u) {
	real_t value = 0.0;
	for (int_t i = 0; i < n; ++i)
		value += v[i]*v[i];
	return 0.5*value;
}

static inline void apply_step_op(real_t ds, real_t dt, complex_t* out, gss_potential_func* potential, grid g, complex_t* wavefunction) {
	for (int_t i = 0; i < g.total_pointcount; ++i) {
		real_t potval = potential(&g.points[g.dimensions*i], g.dimensions, wavefunction[i]);
		real_t opval  = exp(-potval*dt);
		out[i] = wavefunction[i] * opval * ds;
	}
}


























gss_result item_execute(gss_settings settings, gss_potential_func* potential, gss_guess_func* guess) {
	//const real_t Lx = settings.length_x;
	//const real_t Ly = settings.length_y;
	//const real_t dx = settings.length_x / settings.resolution;
	//const real_t dy = settings.length_y / settings.resolution;
	//const real_t ds = (dx*dy);
	//const int_t N = settings.resolution;

	// applying fft -> ifft in fftw scales input by N = n0*n1*...*nk.
	// this is needed to cancel that scaling.
	const real_t ifft_factor = 1.0/(settings.g.total_pointcount);
	const real_t dt = settings.dt;
	real_t ds = 1.0;
	real_t density = 1.0;
	for (int_t i = 0; i < settings.g.dimensions; ++i) {
		ds *= settings.g.deltas[i];
		density *= 1.0/(settings.g.maxs[i] - settings.g.mins[i]);
	}

	gss_result result = {
		.settings = settings,
		.wavefunction = malloc(settings.g.total_pointcount*sizeof(complex_t)),
		.error = 0.0,
		.iterations = 0,
	};

	//real_t* X = result.X;
	//real_t* Y = result.Y;
	//real_t KX[settings.grid.total_pointcount];
	//real_t KY[settings.grid.total_pointcount];
	//FOREACH_ROW(N,N, i,j) {
	//	int_t idx = mat_idx(N,i,j);
	//	X[idx] = -settings.length_x/2 + j*dx;
	//	Y[idx] = -settings.length_y/2 + i*dy;
	//	KX[idx] = (j<N/2) ? 2*M_PI/settings.length_x*j : 2*M_PI/settings.length_x*(-N+j);
	//	KY[idx] = (i<N/2) ? 2*M_PI/settings.length_y*i : 2*M_PI/settings.length_y*(-N+i);
	//}

	grid kgrid = mimic_grid(settings.g);

	//{
	//	int_t indices[settings.g.dimensions];
	//	memset(indices, 0, settings.g.dimensions*sizeof(int_t));
	//	for (int_t i = 0; i < kgrid.total_pointcount; ++i) {
	//		for (int_t j = 0; j < kgrid.dimensions; ++j) {
	//			kgrid.points[kgrid.dimensions*i + j] = 
	//				(indices[j] < kgrid.pointcounts[j]/2) ? 
	//					2*M_PI/(kgrid.maxs[j] - kgrid.mins[j]) * indices[j]
	//				:
	//					2*M_PI/(kgrid.maxs[j] - kgrid.mins[j]) * (-kgrid.pointcounts[j] + indices[j]);
	//		}

	//		int_t baseidx = kgrid.dimensions-1;
	//		indices[baseidx] = (indices[baseidx]+1) % kgrid.pointcounts[baseidx];
	//		for (int_t j = baseidx; j > 0; --j)
	//			if (indices[j] % kgrid.pointcounts[j] == 0)
	//				indices[j-1] = (indices[j-1]+1) % kgrid.pointcounts[j-1];
	//			else
	//				break;
	//	}
	//}

	{
		int_t indices[kgrid.dimensions];
		memset(indices, 0, kgrid.dimensions*sizeof(int_t));
		for (int_t index = 0; index < kgrid.total_pointcount; ++index) {
			int_t prodlen = 1;
			for (int_t n = kgrid.dimensions-1; n >= 0; --n) {
				indices[n] = fmod((index / prodlen), kgrid.pointcounts[n]);
				prodlen *= kgrid.pointcounts[n];
			}

			for (int_t n = 0; n < kgrid.dimensions; ++n) {
				if (indices[n] < kgrid.pointcounts[n]/2) {
					kgrid.points[kgrid.dimensions*index + n] = 2*M_PI / kgrid.lens[n] * indices[n];
				} else {
					kgrid.points[kgrid.dimensions*index + n] = 2*M_PI / kgrid.lens[n] * (-kgrid.pointcounts[n] + indices[n]);
				}
			}
		}
	}

	// Variables used in computation
	//fftw_complex U[settings.grid.total_pointcount];
	fftw_complex old_wavefunction[settings.g.total_pointcount];

	// FFT def.
	complex_t* fft_in   = fftw_malloc(sizeof(complex_t) * settings.g.total_pointcount);
	complex_t* fft_out  = fftw_malloc(sizeof(complex_t) * settings.g.total_pointcount);
	complex_t* ifft_in  = fftw_malloc(sizeof(complex_t) * settings.g.total_pointcount);
	complex_t* ifft_out = fftw_malloc(sizeof(complex_t) * settings.g.total_pointcount);

	fftw_plan fft_plan  = fftw_plan_dft(settings.g.dimensions, settings.g.pointcounts,  fft_in,  fft_out, FFTW_FORWARD,  FFTW_MEASURE);
	fftw_plan ifft_plan = fftw_plan_dft(settings.g.dimensions, settings.g.pointcounts, ifft_in, ifft_out, FFTW_BACKWARD, FFTW_MEASURE);

	// W.f. def.
	{
		//FOREACH_ROW(N,N, i,j) {
		//	int_t idx = mat_idx(N,i,j);
		//	real_t x = X[idx];
		//	real_t y = Y[idx];
		//	result.wavefunction[idx] = guess(x,y);
		//}
		for (int_t i = 0; i < settings.g.total_pointcount; ++i) {
			result.wavefunction[i] = guess(&settings.g.points[settings.g.dimensions*i], settings.g.dimensions);
		}

		real_t sum = 0.0;
		//FOREACH_ROW(N,N, i,j) {
		//	int_t idx = mat_idx(N,i,j);
		//	real_t tmp = cabs(result.wavefunction[idx]);
		//	sum += tmp*tmp;
		//}
		for (int_t i = 0; i < settings.g.total_pointcount; ++i) {
			real_t tmp = cabs(result.wavefunction[i]);
			sum  += tmp*tmp;
		}

		real_t scaling = 1.0/sqrt(sum*ds);
		//FOREACH_ROW(N,N, i,j) {
		//	int_t idx = mat_idx(N,i,j);
		//	result.wavefunction[idx] *= scaling;
		//}
		for (int_t i = 0; i < settings.g.total_pointcount; ++i) {
			result.wavefunction[i] *= scaling;
		}
	}

	apply_step_op(1.0, dt/2.0, fft_in, potential, settings.g, result.wavefunction);
	fftw_execute(fft_plan);

	apply_step_op(ds, dt, ifft_in, VK, kgrid, fft_out);
	fftw_execute(ifft_plan);

	//FOREACH_ROW(N,N, i,j) {
	//	int_t idx = mat_idx(N,i,j);
	//	result.wavefunction[idx] = ifft_factor*ifft_out[idx]/(Lx*Ly);
	//}
	for (int_t i = 0; i < settings.g.total_pointcount; ++i) {
		result.wavefunction[i] = ifft_factor * ifft_out[i] * density;
	}

	for (; result.iterations < settings.max_iterations; ++result.iterations)  {
		memcpy(old_wavefunction, result.wavefunction, settings.g.total_pointcount*sizeof(fftw_complex));

		// Apply operators
		{
			apply_step_op(1.0, dt, fft_in, potential, settings.g, result.wavefunction);
			fftw_execute(fft_plan);

			apply_step_op(ds, dt, ifft_in, VK, kgrid, fft_out);
			fftw_execute(ifft_plan);

			//FOREACH_ROW(N,N, i,j) {
			//	int_t idx = mat_idx(N,i,j);
			//	result.wavefunction[idx] = ifft_factor*ifft_out[idx]/(Lx*Ly);
			//}
			for (int_t i = 0; i < settings.g.total_pointcount; ++i) {
				result.wavefunction[i] = ifft_factor * ifft_out[i] * density;
			}
		}
		
		// Normalize wavefunction
		{
			real_t sum = 0.0;
			//FOREACH_ROW(N,N, i,j) {
			//	int_t idx = mat_idx(N,i,j);
			//	real_t tmp = cabs(result.wavefunction[idx]);
			//	sum += tmp*tmp;
			//}
			for (int_t i = 0; i < settings.g.total_pointcount; ++i) {
				real_t tmp = cabs(result.wavefunction[i]);
				sum += tmp*tmp;
			}

			real_t scaling = 1.0/sqrt(sum*ds);
			//FOREACH_ROW(N,N, i,j) {
			//	int_t idx = mat_idx(N,i,j);
			//	result.wavefunction[idx] *= scaling;
			//}
			for (int_t i = 0; i < settings.g.total_pointcount; ++i) {
				result.wavefunction[i] *= scaling;
			}
		}

		// Calculate error
		{
			real_t sum = 0.0;
			//FOREACH_ROW(N,N, i,j) {
			//	int_t idx = mat_idx(N,i,j);
			//	real_t diff = cabs(result.wavefunction[idx]-old_wavefunction[idx]);
			//	sum += diff*diff;
			//}
			for (int_t i = 0; i < settings.g.total_pointcount; ++i) {
				real_t diff = cabs(result.wavefunction[i] - old_wavefunction[i]);
				sum += diff*diff;
			}
			result.error = sqrt(sum*ds);
		}

		if (settings.measure_every > 0 && result.iterations % settings.measure_every == 0) {
			if (settings.dbgcallback) {
				// Note: old_wavefunction has already been used at this point, we can
				// thus use it as a temporary buffer.
				apply_step_op(1.0, dt/2.0, old_wavefunction, potential, settings.g, result.wavefunction);
				settings.dbgcallback(settings.g, old_wavefunction);
			}
		}

		if (result.error < settings.error_tol) {
			break;
		}
	}

	apply_step_op(1.0, dt/2.0, result.wavefunction, potential, settings.g, result.wavefunction);

	// Cleanup
	free_grid(kgrid);
	fftw_destroy_plan(fft_plan);
	fftw_destroy_plan(ifft_plan);
	fftw_free(fft_in);
	fftw_free(fft_out);
	fftw_free(ifft_in);
	fftw_free(ifft_out);

	return result;
}
