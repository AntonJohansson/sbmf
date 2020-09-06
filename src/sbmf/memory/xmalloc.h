#pragma once

#include <sbmf/types.h>
#include <sbmf/debug/log.h>

#include <stdlib.h>

static inline void* _xmalloc_impl(u64 size_in_bytes, const u32 linenumber,
																	const char file[], const char func[]) {
	/*
	log_info("[%s:%d in %s()] allocating %ld bytes via malloc.", file, linenumber,
					 func, size_in_bytes);
	*/
	void* mem = malloc(size_in_bytes);
	if(!mem) {
		log_error("malloc(...) could not allocate %lu bytes.", size_in_bytes);
	}

	return mem;
}

#define xmalloc(size_in_bytes) \
	_xmalloc_impl(size_in_bytes, __LINE__, __FILE__, __func__)
