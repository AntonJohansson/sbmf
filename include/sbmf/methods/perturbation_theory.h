#pragma once

#include <sbmf/methods/nlse_solver.h>

struct pt_result {
	f64 E0, E1, E2;
};

struct pt_result rayleigh_schroedinger_pt(struct nlse_result res, f64* g0, u32* particle_count);
