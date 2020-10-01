#include "find_groundstate.h"

#include <fftw3.h>

#include <string.h> /* memcpy */
#include <assert.h>

static f64 VK(f64* v, i32 n, c64 u) {
	SBMF_UNUSED(u);
	f64 value = 0.0;
	for (i32 i = 0; i < n; ++i)
		value += v[i]*v[i];
	return 0.5*value;
}

static void apply_step_op(f64 ds, f64 dt, c64* out, gss_potential_func* potential, struct grid g, c64* wavefunction) {
	for (i32 i = 0; i < g.total_pointcount; ++i) {
		f64 potval = potential(&g.points[g.dimensions*i], g.dimensions, wavefunction[i]);
		f64 opval  = exp(-potval*dt);
		out[i] = wavefunction[i] * opval * ds;
	}
}

struct gss_result item(struct item_settings settings, gss_potential_func* potential, gss_guess_func* guess) {
	/*
	 * applying fft -> ifft in fftw scales input by N = n0*n1*...*nk.
	 * this is needed to cancel that scaling.
	 */
	const f64 ifft_factor = 1.0/(settings.g.total_pointcount);
	const f64 dt = settings.dt;
	f64 ds = 1.0;
	f64 density = 1.0;
	for (i32 i = 0; i < settings.g.dimensions; ++i) {
		ds *= settings.g.deltas[i];
		density *= 1.0/(settings.g.maxs[i] - settings.g.mins[i]);
	}

	struct gss_result result = {
		.wavefunction = (c64*) sbmf_stack_push(settings.g.total_pointcount*sizeof(c64)),
		.error = 0.0,
		.iterations = 0,
	};

	struct grid kgrid = mimic_grid(settings.g);

	/* Momentum grid */
	{
		i32 indices[kgrid.dimensions];
		memset(indices, 0, kgrid.dimensions*sizeof(i32));
		for (i32 index = 0; index < kgrid.total_pointcount; ++index) {
			i32 prodlen = 1;
			for (i32 n = kgrid.dimensions-1; n >= 0; --n) {
				indices[n] = fmod((index / prodlen), kgrid.pointcounts[n]);
				prodlen *= kgrid.pointcounts[n];
			}

			for (i32 n = 0; n < kgrid.dimensions; ++n) {
				if (indices[n] < kgrid.pointcounts[n]/2) {
					kgrid.points[kgrid.dimensions*index + n] = 2*M_PI / kgrid.lens[n] * indices[n];
				} else {
					kgrid.points[kgrid.dimensions*index + n] = 2*M_PI / kgrid.lens[n] * (-kgrid.pointcounts[n] + indices[n]);
				}
			}
		}
	}

	/* Variables used in computation */
	fftw_complex old_wavefunction[settings.g.total_pointcount];

	/* FFT def. */
	c64* fft_in   = fftw_malloc(sizeof(c64) * settings.g.total_pointcount);
	c64* fft_out  = fftw_malloc(sizeof(c64) * settings.g.total_pointcount);
	c64* ifft_in  = fftw_malloc(sizeof(c64) * settings.g.total_pointcount);
	c64* ifft_out = fftw_malloc(sizeof(c64) * settings.g.total_pointcount);

	fftw_plan fft_plan  = fftw_plan_dft(settings.g.dimensions, settings.g.pointcounts,  fft_in,  fft_out, FFTW_FORWARD,  FFTW_MEASURE);
	fftw_plan ifft_plan = fftw_plan_dft(settings.g.dimensions, settings.g.pointcounts, ifft_in, ifft_out, FFTW_BACKWARD, FFTW_MEASURE);

	/* Spatial grid */
	{
		for (i32 i = 0; i < settings.g.total_pointcount; ++i) {
			result.wavefunction[i] = guess(&settings.g.points[settings.g.dimensions*i], settings.g.dimensions);
		}

		f64 sum = 0.0;
		for (i32 i = 0; i < settings.g.total_pointcount; ++i) {
			f64 tmp = cabs(result.wavefunction[i]);
			sum  += tmp*tmp;
		}

		f64 scaling = 1.0/sqrt(sum*ds);
		for (i32 i = 0; i < settings.g.total_pointcount; ++i) {
			result.wavefunction[i] *= scaling;
		}
	}

	apply_step_op(1.0, dt/2.0, fft_in, potential, settings.g, result.wavefunction);
	fftw_execute(fft_plan);

	apply_step_op(ds, dt, ifft_in, VK, kgrid, fft_out);
	fftw_execute(ifft_plan);

	for (i32 i = 0; i < settings.g.total_pointcount; ++i) {
		result.wavefunction[i] = ifft_factor * ifft_out[i] * density;
	}

	for (; result.iterations < settings.max_iterations; ++result.iterations)  {
		memcpy(old_wavefunction, result.wavefunction, settings.g.total_pointcount*sizeof(fftw_complex));

		/* Apply operators */
		{
			apply_step_op(1.0, dt, fft_in, potential, settings.g, result.wavefunction);
			fftw_execute(fft_plan);

			apply_step_op(ds, dt, ifft_in, VK, kgrid, fft_out);
			fftw_execute(ifft_plan);

			for (i32 i = 0; i < settings.g.total_pointcount; ++i) {
				result.wavefunction[i] = ifft_factor * ifft_out[i] * density;
			}
		}

		/* Normalize wavefunction */
		{
			f64 sum = 0.0;
			for (i32 i = 0; i < settings.g.total_pointcount; ++i) {
				f64 tmp = cabs(result.wavefunction[i]);
				sum += tmp*tmp;
			}

			f64 scaling = 1.0/sqrt(sum*ds);
			for (i32 i = 0; i < settings.g.total_pointcount; ++i) {
				result.wavefunction[i] *= scaling;
			}
		}

		/* Calculate error */
		{
			f64 sum = 0.0;
			for (i32 i = 0; i < settings.g.total_pointcount; ++i) {
				f64 diff = cabs(result.wavefunction[i] - old_wavefunction[i]);
				sum += diff*diff;
			}
			result.error = sqrt(sum*ds);
		}

		if (settings.measure_every > 0 && result.iterations % settings.measure_every == 0) {
			if (settings.dbgcallback) {
				/*
				 * Note: old_wavefunction has already been used at this point, we can
				 * thus use it as a temporary buffer.
				 */
				apply_step_op(1.0, dt/2.0, old_wavefunction, potential, settings.g, result.wavefunction);
				settings.dbgcallback(settings, old_wavefunction);
			}
		}

		if (result.error < settings.error_tol) {
			break;
		}
	}

	apply_step_op(1.0, dt/2.0, result.wavefunction, potential, settings.g, result.wavefunction);

	/*
	 * Cleanup
	 * This is not handledefasdfasf free_grid(kgrid);
	 */
	fftw_destroy_plan(fft_plan);
	fftw_destroy_plan(ifft_plan);
	fftw_free(fft_in);
	fftw_free(fft_out);
	fftw_free(ifft_in);
	fftw_free(ifft_out);

	return result;
}
