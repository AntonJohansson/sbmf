static inline void* xmalloc(u64 size_in_bytes) {
	void* mem = malloc(size_in_bytes);
	if(!mem) {
		sbmf_log_panic("malloc(...) could not allocate %lu bytes.", size_in_bytes);
	}

	return mem;
}

/*
 * Very basic, non-expanding stack-allocator
 */

struct stack_allocator {
	u32 top;
	u32 size;
	u8* memory;
};

static struct stack_allocator* sa_make(u32 size_in_bytes) {
	void* mem = xmalloc(sizeof(struct stack_allocator) + size_in_bytes);

	struct stack_allocator* sa = (struct stack_allocator*) mem;
	sa->top = 0;
	sa->size = size_in_bytes;
	sa->memory = (u8*) (sa+1);

	return sa;
}

static void sa_destroy(struct stack_allocator* sa) {
	free(sa);
}

static u8* sa_push(struct stack_allocator* sa, u32 size_in_bytes) {
	if (sa->top + size_in_bytes > sa->size) {
		sbmf_log_panic("stack allocator %p out of memory (%u/%u bytes)\n\ttrying to allocate %u bytes\n", sa, sa->top, sa->size, size_in_bytes);
	}

	u32 alignment_padding = 8 - (sa->top % 8);
	sa->top += alignment_padding;

	u8* ptr = sa->memory + sa->top;
	sa->top += size_in_bytes;

	return ptr;
}
