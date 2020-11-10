#pragma once

#include <stdlib.h>

static inline void* xmalloc(u64 size_in_bytes) {
	void* mem = malloc(size_in_bytes);
	if(!mem) {
		sbmf_log_panic("malloc(...) could not allocate %lu bytes.", size_in_bytes);
	}

	return mem;
}
