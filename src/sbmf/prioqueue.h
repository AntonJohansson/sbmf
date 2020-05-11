#pragma once

#include "common/common.h"
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <string.h> // memcpy, memmove

// Returns true if a comes before b
typedef bool cmpfunc(void* a, void* b);

typedef struct {
	i32 top;
	u32 size;

	u32 max_size;
	u32 item_size;

	u8* mem;

	cmpfunc* cmp;
} prioqueue;

static inline prioqueue* prioqueue_new(u32 items, u32 item_size, cmpfunc* cmp) {
	void* mem = malloc(items*item_size + sizeof(prioqueue));
	assert(mem);

	prioqueue* pq = (prioqueue*) mem;

	pq->top = 0;
	pq->size = 0;

	pq->max_size = items;
	pq->item_size = item_size;

	pq->mem = (u8*) (pq + 1);

	pq->cmp = cmp;

	return pq;
}

static inline void prioqueue_free(prioqueue* pq) {
	free(pq);
}

static inline void* prioqueue_get(prioqueue* pq, i32 index) {
	return &pq->mem[((pq->top + index) % pq->max_size) * pq->item_size];
}

static inline void prioqueue_push(prioqueue* pq, void* data) {
	assert(pq->size < pq->max_size);

	i32 insert_index = 0;
	for (; insert_index < pq->size; ++insert_index) {
		if (pq->cmp(data, prioqueue_get(pq, insert_index)))
			break;
	}

	i32 wrapped_index = (pq->top + insert_index) % pq->max_size;

	if (insert_index < pq->size)  {
				// 0 1 2 3 4 5 6 7 8 9
		// 	   b   t     i
		// shift 012 -> 123
		// shift 9 -> 0
		// shift 78 -> 89

		i32 end_index = (pq->top + pq->size) % pq->max_size;

		if (wrapped_index < end_index) {
			// 0 1 2 3 4 5 6 7 8 9
			// 	   t   i
			// shift 456 -> 567

			i32 next_index = (wrapped_index + 1) % pq->max_size;
			memmove(&pq->mem[next_index*pq->item_size], &pq->mem[wrapped_index*pq->item_size], (pq->size - insert_index)*pq->item_size);
		} else {
			memmove(&pq->mem[pq->item_size], &pq->mem[0], (end_index)*pq->item_size);

			memmove(&pq->mem[0], &pq->mem[(pq->max_size-1)*pq->item_size], pq->item_size);

			i32 next_index = (wrapped_index + 1) % pq->max_size;
			memmove(&pq->mem[next_index*pq->item_size], &pq->mem[wrapped_index*pq->item_size], (pq->size - insert_index-1-end_index)*pq->item_size);
		}
	}

	memcpy(&pq->mem[wrapped_index*pq->item_size], data, pq->item_size);
	pq->size++;
}

static inline u8* prioqueue_top(prioqueue* pq) {
	return &pq->mem[pq->top * pq->item_size];
}

static inline bool prioqueue_pop(prioqueue* pq) {
	if (pq->size > 0) {
		pq->top = (pq->top + 1) % pq->max_size;
		pq->size--;

		return (pq->size > 0);
	}

	return false;
}

