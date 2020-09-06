#pragma once

#include "types.h"

void sbmf_init();
void sbmf_shutdown();

u8* sbmf_stack_push(u32 size_in_bytes);
u32 sbmf_stack_marker();
void sbmf_stack_free_to_marker(u32 marker);

#define xmalloc(size_in_bytes) \
	_xmalloc_impl(size_in_bytes, __LINE__, __FILE__, __func__)
