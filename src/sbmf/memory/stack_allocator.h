#pragma once

#include <sbmf/common/common.h>
//#include <sbmf/memory/memory.h>
#include <sbmf/memory/xmalloc.h>
#include <sbmf/common/log.h>

#define SA_PUSH(sa, type) \
	(type*)sa_push(sa, sizeof(type))

typedef struct stack_allocator stack_allocator;
struct stack_allocator {
	u32 top;
	u32 size;
	u8* memory;
};

static inline stack_allocator* sa_make(u32 size_in_bytes) {
	void* mem = xmalloc(sizeof(stack_allocator) + size_in_bytes);

	stack_allocator* sa = (stack_allocator*) mem;
	sa->top = 0;
	sa->size = size_in_bytes;
	sa->memory = (u8*) (sa+1);

	return sa;
}

static inline void sa_destroy(stack_allocator* sa) {
	free(sa);
}

static inline u8* sa_push(stack_allocator* sa, u32 size_in_bytes) {
	log_info("stack allocator %p (%u/%u bytes)\n\ttrying to allocate %u bytes\n", sa, sa->top, sa->size, size_in_bytes);
	if (sa->top + size_in_bytes > sa->size) {
		log_error("stack allocator %p out of memory (%u/%u bytes)\n\ttrying to allocate %u bytes\n", sa, sa->top, sa->size, size_in_bytes);
		return 0;
	}

	u8* ptr = sa->memory + sa->top;
	sa->top += size_in_bytes;

	return ptr;
}
