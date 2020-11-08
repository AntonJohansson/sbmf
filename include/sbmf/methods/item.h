#pragma once

#include <sbmf/types.h>
#include <sbmf/math/grid.h>

struct item_settings;

typedef void item_debug_callback(struct item_settings, c64*);
typedef f64 gss_potential_func(f64* v, i32 n, c64 u);
typedef c64  gss_guess_func(f64* v, i32 n);

struct item_settings {
	struct grid g;

	u32 max_iterations;
	f64 error_tol;

	u32 measure_every;
	item_debug_callback* dbgcallback;

	f64 dt;
};

struct gss_result {
	c64* wavefunction;

	f64 error;
	u32 iterations;
};

struct gss_result item(struct item_settings settings, gss_potential_func* potential, gss_guess_func* guess);
