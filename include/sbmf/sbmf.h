#pragma once

#include "types.h"

#define SBMF_UNUSED(x) \
	(void)x

void sbmf_init();
void sbmf_shutdown();

u8* sbmf_stack_push_impl(u32 size_in_bytes, const u32 linenumber, const char file[], const char func[]);
#define sbmf_stack_push(size_in_bytes) \
	sbmf_stack_push_impl(size_in_bytes, __LINE__, __FILE__, __func__)

u32 sbmf_stack_marker();
void sbmf_stack_free_to_marker(u32 marker);
