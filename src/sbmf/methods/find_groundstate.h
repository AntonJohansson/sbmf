#pragma once

#include <sbmf/types.h>
#include <sbmf/math/grid.h>

struct item_settings;
struct scim_settings;
struct gss_result;

typedef f64 gss_potential_func(f64* v, i32 n, c64 u);

/* Does not handle 2D/3D/.. */
typedef void gss_potential_vec_func(f64* out, f64* in_x, c64* in_u, u32 len);

typedef c64  gss_guess_func(f64* v, i32 n);
typedef void gss_guess_vec_func(c64* out, f64* in_x, u32 len);

typedef void item_debug_callback(struct item_settings, c64*);
typedef void scim_debug_callback(struct scim_settings, c64*);

struct item_settings {
	struct grid g;

	u32 max_iterations;
	f64 error_tol;

	u32 measure_every;
	item_debug_callback* dbgcallback;

	f64 dt;
};

struct scim_settings {
	u32 num_basis_functions;

	u32 max_iterations;
	f64 error_tol;

	u32 measure_every;
	scim_debug_callback* dbgcallback;
};

struct gss_result {
	c64* wavefunction;

	f64 error;
	u32 iterations;
};


struct gss_result item(struct item_settings settings, gss_potential_func* potential, gss_guess_func* guess);
struct gss_result scim(struct scim_settings settings, gss_potential_vec_func* potential, gss_guess_vec_func* guess);
