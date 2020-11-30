#pragma once

#include <sbmf/methods/nlse_solver.h>

struct bestmf_result {
	f64 energy;
	u32 coeff_count;
	u32 comp_count;
	f64* coeff;
	f64 n1;
	f64 n2;
};

struct bestmf_result best_meanfield(struct nlse_settings settings, const u32 particle_count, f64 g0, struct nlse_guess* guesses);

struct bestmf_2comp_result {
	f64 energy;
	u32 coeff_count;
	u32 comp_count;
	f64* coeff;
	f64 n1;
	f64 n2;
	f64 n3;
	f64 n4;
};
struct bestmf_2comp_result best_meanfield_2comp(struct nlse_settings settings,
		const u32 particle_count,
		f64 g0[static 2*2],
		struct nlse_guess* guesses);
