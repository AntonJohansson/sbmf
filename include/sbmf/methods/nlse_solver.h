#pragma once

#include <sbmf/types.h>
#include <sbmf/math/matrix.h>
#include <sbmf/methods/quadgk_vec.h>
#include <sbmf/math/basis.h>

/*
 * Solves a general, non-linear system of
 * Schrödinger equations by the means of
 * iteration until self-consistency in
 * a specified basis.
 */

typedef void nlse_operator_func(const u32 len, f64 out[static len],
		f64 in_x[static len], const u32 component_count,
		f64 in_u[static len*component_count],
		void* userdata);
typedef void nlse_callback(c64* a, c64* b, u32 len);

struct nlse_settings;
struct nlse_result;
typedef void nlse_debug_callback(struct nlse_settings, struct nlse_result);

struct nlse_settings {
	u32 max_iterations;
	f64 error_tol;

	/* Separating out the spatial potentiential
	 * allows for optimizations */
	nlse_operator_func* spatial_pot_perturbation;

	u32 measure_every;
	nlse_debug_callback* debug_callback;

	nlse_callback* post_normalize_callback;

	/* Choice of GK rule used for numerical integration internally */
	struct gk_data gk;

	/* Choice of basis to solve to problem in */
	u32 num_basis_funcs;
	struct basis basis;

	/* Everything below the zero_threshold is considered
	 * 0 in the hamiltonian. */
	f64 zero_threshold;
};

/* Initial Guess */
typedef void nlse_coeff_guess_func(f64* out, u32 len, u32 component);
typedef void nlse_spatial_guess_func(f64* out, f64* in, u32 len, void* data);

struct nlse_guess {
	enum {
		DEFAULT_GUESS 	= 0,
		SPATIAL_GUESS 	= 1,
		COEFF_GUESS   	= 2
	} type;
	union {
		nlse_coeff_guess_func* 	 coeff_guess;
		nlse_spatial_guess_func* spatial_guess;
	} data;
};

struct nlse_component {
	struct nlse_guess guess;
	nlse_operator_func* op;
	void* userdata;
};

struct nlse_result {
	u32 iterations;
	u32 component_count;
	u32 coeff_count;
	f64* coeff;
	f64* error;
	f64* energy;
	struct hermitian_bandmat* hamiltonian;
};

struct nlse_result nlse_solver(struct nlse_settings settings, const u32 component_count, struct nlse_component components[static component_count]);
