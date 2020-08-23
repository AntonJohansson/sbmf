#pragma once

#include <sbmf/types.h>
#include <sbmf/math/grid.h>

typedef struct gss_settings gss_settings;
typedef struct gss_result gss_result;

typedef f64 gss_potential_func(f64* v, i32 n, c64 u);
typedef c64  gss_guess_func(f64* v, i32 n);
typedef void gss_debug_callback(struct grid g, c64* wf);

typedef struct {
	u32 n[2];
	c64* coeffs;
	u32 coeff_count;
	gss_potential_func* pot;
} hob_integrand_params;

struct gss_settings {
	struct grid g;

	u32 max_iterations;
	f64 error_tol;

	u32 measure_every;
	gss_debug_callback* dbgcallback;

	f64 dt;
};

struct gss_result {
	gss_settings settings;

	c64* wavefunction;

	f64 error;
	u32 iterations;
};

void gss_free_result(gss_result res);
gss_result item_execute(gss_settings settings, gss_potential_func* potential, gss_guess_func* guess);
gss_result hob(gss_settings settings, gss_potential_func* potential, gss_guess_func* guess);

void hob_perf_test();
