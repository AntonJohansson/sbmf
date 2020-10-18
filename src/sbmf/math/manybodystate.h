#pragma once

#include <sbmf/types.h>

/*
 * |0> = |N,0,...>
 * |mn> = |N-2,0,...,0,1,0,...,0,1,0,...>
 * 					   ^-_ m:th  ^-_ n:th
 * |> = c0|...> + c1|...> + ...
 * <x|> = c0<x|...> + c1<x|...> + ...
 * <x|abc...> = phi_a(x)phi_b(x)phi_c(x)...
 */
struct manybodystate {
	u32 index_count;
	struct {
		u32 particle_count;
		u32 index_in_basis;
	}* single_states;
};

/* How do we calculate the energy of a manybody state? */
