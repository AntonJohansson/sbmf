#pragma once

#include <sbmf/common/common.h>
#include <sbmf/common/grid.h>

struct gss_settings;
typedef struct gss_settings gss_settings;
struct gss_result;
typedef struct gss_result gss_result;

typedef real_t gss_potential_func(real_t* v, int_t n, complex_t u);
typedef complex_t  gss_guess_func(real_t* v, int_t n);
typedef void gss_debug_callback(grid g, complex_t* wf);

struct gss_settings {
	grid g;

	int_t max_iterations;
	real_t error_tol;

	int_t measure_every;
	gss_debug_callback* dbgcallback;

	// Specific for aitem
	real_t c;
	real_t A;
	real_t dt;
};

struct gss_result {
	gss_settings settings;

	complex_t* wavefunction;

	real_t error;
	int_t iterations;
};

extern void gss_free_result(gss_result res);

extern gss_result aitem_execute(gss_settings settings,
																gss_potential_func* potential, gss_guess_func* guess);

extern gss_result niter_execute(gss_settings settings,
																gss_potential_func* potential, gss_guess_func* guess);
extern gss_result item_execute(gss_settings settings,
																gss_potential_func* potential, gss_guess_func* guess);
