#pragma once

#include <sbmf/types.h>

void find_best_meanfield_occupations(const u32 particle_count,
		const u32 coeff_count,
		c64 orbital_1_coeffs[static coeff_count],
		c64 orbital_2_coeffs[static coeff_count]);
