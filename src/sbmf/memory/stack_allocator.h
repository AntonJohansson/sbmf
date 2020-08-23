#pragma once

#include <sbmf/types.h>
#include <sbmf/memory/xmalloc.h>
#include <assert.h>

#ifndef SA_PRINT_DEBUG_INFO
#define SA_PRINT_DEBUG_INFO 0
#endif

#define SA_PUSH(sa, type) \
	(type*)sa_push(sa, sizeof(type))

struct stack_allocator {
	u32 top;
	u32 size;
	u8* memory;
};

static inline struct stack_allocator* sa_make(u32 size_in_bytes) {
	void* mem = xmalloc(sizeof(struct stack_allocator) + size_in_bytes);

	struct stack_allocator* sa = (struct stack_allocator*) mem;
	sa->top = 0;
	sa->size = size_in_bytes;
	sa->memory = (u8*) (sa+1);

	return sa;
}

static inline void sa_destroy(struct stack_allocator* sa) {
	free(sa);
}

static inline u8* sa_push(struct stack_allocator* sa, u32 size_in_bytes) {
#if SA_PRINT_DEBUG_INFO
	log_info("stack allocator %p (%u/%u bytes)\n\ttrying to allocate %u bytes\n", sa, sa->top, sa->size, size_in_bytes);
#endif

	if (sa->top + size_in_bytes > sa->size) {
		log_panic("stack allocator %p out of memory (%u/%u bytes)\n\ttrying to allocate %u bytes\n", sa, sa->top, sa->size, size_in_bytes);
	}

	u32 alignment_padding = 8 - (sa->top % 8);
	sa->top += alignment_padding;

	u8* ptr = sa->memory + sa->top;
	sa->top += size_in_bytes;

	return ptr;
}
