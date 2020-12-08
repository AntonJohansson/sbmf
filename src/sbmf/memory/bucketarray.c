#pragma once

#include <sbmf/sbmf.h>
#include <string.h>

struct bucket_header {
	struct bucket_header* next;
	u8* mem;
};

struct barray {
	u32 items_per_bucket;
	u32 item_size;

	u32 bucket_count;

	struct bucket_header* base;
};

static inline struct bucket_header* barray_new_bucket(struct barray* ba) {
	struct bucket_header* hdr = (struct bucket_header*)sbmf_stack_push(sizeof(struct bucket_header) + ba->items_per_bucket*ba->item_size);

	hdr->next = 0;
	hdr->mem = (u8*) (hdr + 1);

	ba->bucket_count++;

	return hdr;
}

static inline struct barray* barray_new(u32 items_per_bucket, u32 item_size) {
	struct barray* ba = (struct barray*) sbmf_stack_push(sizeof(struct barray));
	ba->items_per_bucket = items_per_bucket;
	ba->item_size = item_size;
	ba->bucket_count = 0;
	ba->base = barray_new_bucket(ba);

	return ba;
}

static inline void* barray_get(struct barray* ba, u32 index) {
	u32 bucket = index / ba->items_per_bucket;
	u32 index_in_bucket = index % ba->items_per_bucket;

	/* Appends bucket until it finds the desired bucket. */
	struct bucket_header* ptr = ba->base;
	while (bucket--) {
		if (!ptr->next)
			ptr->next = barray_new_bucket(ba);

		ptr = ptr->next;
	}

	return &ptr->mem[index_in_bucket * ba->item_size];
}

static inline void barray_set(struct barray* ba, u32 index, const void* item) {
	void* ptr = barray_get(ba, index);
	memcpy(ptr, item, ba->item_size);
}
