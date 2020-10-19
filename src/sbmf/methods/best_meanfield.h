#pragma once

#include <sbmf/types.h>

void find_best_meanfield_occupations(const u32 particle_count,
									 const f64 g,
                                     const u32 coeff_count,
                                     c64 orbital_1_coeffs[static coeff_count],
                                     c64 orbital_2_coeffs[static coeff_count],
                                     const f64 energy_1,
                                     const f64 energy_2);

/* Assumes our manybody states consists of two orbitals phi1,phi2 and occupation numbers n1,n2 */
void best_meanfield_energy(const u32 particle_count,
                           const u32 coeff_count,
                           c64 orbital_1_coeffs[static coeff_count],
                           c64 orbital_2_coeffs[static coeff_count],
						   const f64 energy_1,
						   const f64 energy_2,
						   const u32 occupation_1,
						   const u32 occupation_2);
