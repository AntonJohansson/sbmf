#pragma once

#include <sbmf/types.h>
#include <sbmf/math/grid.h>
#include <sbmf/math/matrix.h>

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

typedef void gp2c_operator_func(f64* out, f64* in_x, c64* in_a, c64* in_b, u32 len);
typedef void gp2c_guess_func(c64* out, u32 len);
typedef void gp2c_callback(c64* a, c64* b, u32 len);

struct gp2c_settings {
	u32 num_basis_functions;

	u32 max_iterations;
	f64 error_tol;

	u32 measure_every;
	scim_debug_callback* dbgcallback;

	gp2c_guess_func* guess_a;
	gp2c_guess_func* guess_b;

	gp2c_callback* post_normalize_callback;
};

struct gp2c_result {
	c64* coeff_a;
	c64* coeff_b;

	f64 error_a;
	f64 error_b;
	u32 iterations;

	f64 energy_a;
	f64 energy_b;

	struct complex_hermitian_bandmat hamiltonian_a;
	struct complex_hermitian_bandmat hamiltonian_b;
};

struct gp2c_result gp2c(struct gp2c_settings, gp2c_operator_func* op_a, gp2c_operator_func* op_b);
