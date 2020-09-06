#pragma once

#include <sbmf/types.h>
#include <sbmf/math/grid.h>

struct gss_settings;
struct gss_result;

typedef f64 gss_potential_func(f64* v, i32 n, c64 u);
typedef c64  gss_guess_func(f64* v, i32 n);
typedef void gss_debug_callback(struct gss_settings, c64*);

struct gss_settings {
	struct grid g;

	u32 max_iterations;
	f64 error_tol;

	u32 measure_every;
	gss_debug_callback* dbgcallback;

	f64 dt;
};

struct gss_result {
	struct gss_settings settings;

	c64* wavefunction;

	f64 error;
	u32 iterations;
};

struct gss_result item(struct gss_settings settings, gss_potential_func* potential, gss_guess_func* guess);
struct gss_result scim(struct gss_settings settings, gss_potential_func* potential, gss_guess_func* guess);
