#pragma once

#include <sbmf/types.h>
#include <sbmf/math/grid.h>
#include <sbmf/math/matrix.h>
#include <sbmf/methods/quadgk_vec.h>

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

struct gss_result ho_scim(struct scim_settings settings, gss_potential_vec_func* potential, gss_guess_vec_func* guess);

/* What would a two component GP-system look like?
 * 		Call components a,b
 *			T(a) + V(a) + g_a|a|^2 + g_ab|b|^2
 *			T(b) + V(b) + g_b|b|^2 + g_ba|a|^2
 *		where only really g_a,g_b,g_ab,g_ba are variables
 *		needed to be passed to the solver. All other things can
 *		be inlined for performance.
 */

typedef void gp2c_operator_func(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                c64 in_u[static len*component_count]);
typedef void gp2c_guess_func(c64* out, u32 len);
typedef void gp2c_callback(c64* a, c64* b, u32 len);

struct gp2c_settings;
struct gp2c_result;
typedef void gp2c_debug_callback(struct gp2c_settings, struct gp2c_result);

struct gp2c_settings {
	u32 num_basis_functions;

	u32 max_iterations;
	f64 error_tol;

	gp2c_operator_func* ho_potential_perturbation;

	u32 measure_every;
	gp2c_debug_callback* dbgcallback;

	gp2c_callback* post_normalize_callback;

	/* Choice of GK rule used for numerical integration internally */
	struct gk_data gk;
};

struct gp2c_component {
	gp2c_guess_func* guess;
	gp2c_operator_func* op;
};

struct gp2c_result {
	u32 iterations;
	u32 component_count;
	u32 coeff_count;
	c64* coeff;
	f64* error;
	f64* energy;
	struct complex_hermitian_bandmat* hamiltonian;
};

struct gp2c_result gp2c(struct gp2c_settings settings, const u32 component_count, struct gp2c_component components[static component_count]);
