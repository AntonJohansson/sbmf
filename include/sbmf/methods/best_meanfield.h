#pragma once

#include <sbmf/methods/nlse_solver.h>

struct nlse_result best_meanfield(struct nlse_settings settings, const u32 particle_count, f64 g0);
