#pragma once

#include <sbmf/types.h>
#include <sbmf/math/matrix.h>
#include <sbmf/methods/quadgk_vec.h>
#include <sbmf/math/basis.h>

typedef void gp2c_operator_func(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count],
								void* userdata);
typedef void gp2c_guess_func(f64* out, u32 len);
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

	/* Choice of basis to solve to problem in */
	struct basis basis;
};

struct gp2c_component {
	gp2c_guess_func* guess;
	gp2c_operator_func* op;
	void* userdata;
};

struct gp2c_result {
	u32 iterations;
	u32 component_count;
	u32 coeff_count;
	f64* coeff;
	f64* error;
	f64* energy;
	struct hermitian_bandmat* hamiltonian;
};

struct gp2c_result gp2c(struct gp2c_settings settings, const u32 component_count, struct gp2c_component components[static component_count]);
