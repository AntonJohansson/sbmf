#pragma once

#include <sbmf/methods/nlse_solver.h>

struct gp_settings {
	nlse_operator_func* pot;

	u32 num_basis_funcs;
	struct basis basis;

	u32 component_count;
	u32* occupations; 				/* [static component_count] */
	struct nlse_guess* guesses; 	/* [static component_count] */
	f64* g0; 						/* [static component_count*component_count] */

	nlse_debug_callback* debug_callback;
	u32 measure_every;

	f64 zero_threshold;
};

struct nlse_result grosspitaevskii(struct gp_settings settings);

f64 full_energy(struct gp_settings settings, struct nlse_result res);
