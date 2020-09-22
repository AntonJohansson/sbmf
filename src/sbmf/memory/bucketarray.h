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

	struct bucket_header* base;
};

inline struct bucket_header* barray_add_bucket(struct barray* ba) {
	struct bucket_header* hdr = (struct bucket_header*)sbmf_stack_push(sizeof(struct bucket_header) + ba->items_per_bucket*ba->item_size);

	hdr->next = 0;
	hdr->mem = (u8*) (hdr + 1);

	return hdr;
}

inline struct barray* barray_new(u32 items_per_bucket, u32 item_size) {
	struct barray* ba = (struct barray*) sbmf_stack_push(sizeof(struct barray));
	ba->items_per_bucket = items_per_bucket;
	ba->item_size = item_size;
	ba->base = barray_new_bucket(ba);

	return ba;
}

inline void* barray_get(struct barray* ba, u32 index) {
	u32 bucket = index / ba->items_per_bucket;
	u32 index_in_bucket = index % ba->items_per_bucket;

	/* Appends bucket until it finds the desired bucket. */
	struct barray* ptr = ba->base;
	while (bucket--) {
		if (!ptr->next)
			ptr->next = barray_new_bucket(ba);

		ptr = ptr->next;
	}

	return &ptr->mem[index_in_bucket * ba->item_size];
}

inline void barray_set(struct barray* ba, u32 index, const void* item) {
	void* ptr = barray_get(ba, index);
	memcpy(ptr, item, ba->item_size);
}
