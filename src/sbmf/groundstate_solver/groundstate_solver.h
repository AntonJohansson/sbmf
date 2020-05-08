#pragma once

#include <sbmf/common/common.h>
#include <sbmf/common/grid.h>

struct gss_settings;
typedef struct gss_settings gss_settings;
struct gss_result;
typedef struct gss_result gss_result;

typedef f64 gss_potential_func(f64* v, i32 n, c64 u);
typedef c64  gss_guess_func(f64* v, i32 n);
typedef void gss_debug_callback(grid g, c64* wf);

struct gss_settings {
	grid g;

	i32 max_iterations;
	f64 error_tol;

	i32 measure_every;
	gss_debug_callback* dbgcallback;

	f64 dt;
};

struct gss_result {
	gss_settings settings;

	c64* wavefunction;

	f64 error;
	i32 iterations;
};

extern void gss_free_result(gss_result res);

extern gss_result niter_execute(gss_settings settings,
																gss_potential_func* potential, gss_guess_func* guess);
extern gss_result item_execute(gss_settings settings,
																gss_potential_func* potential, gss_guess_func* guess);
