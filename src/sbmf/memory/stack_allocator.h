#pragma once

#include <sbmf/common/common.h>
#include <sbmf/memory/memory.h>
#include <assert.h>

#define SA_PUSH(sa, type) \
	(type*)sa_push(sa, sizeof(type))

typedef struct {
	u32 top;
	u32 size;
	u8* memory;
} stack_allocator;

static inline stack_allocator* sa_new(u32 size_in_bytes) {
	void* mem = malloc(sizeof(stack_allocator) + size_in_bytes);
	assert(mem);

	stack_allocator* sa = (stack_allocator*) mem;
	sa->top = 0;
	sa->size = size_in_bytes;
	sa->memory = (u8*) (sa+1);

	return sa;
}

static inline stack_allocator sa_free(stack_allocator* sa) {
	free(sa);
}

static inline u8* sa_push(stack_allocator* sa, u32 size_in_bytes) {
	assert(sa->top + size_in_bytes <= sa->size);

	u8* ptr = sa->memory + sa->top;
	sa->top += size_in_bytes;

	return ptr;
}
