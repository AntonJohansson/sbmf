#pragma once

#include <sbmf/sbmf.h>
#include <string.h> /* memcpy, memove */
#include "bucketarray.h"

#define PQ_PARENT(x) (((i32)x-1)/2)
#define PQ_LEFT(x)	 (2*(i32)x+1)
#define PQ_RIGHT(x)  (2*(i32)x+2)
/*
 * | 0 |  1 |  2 |  3 |  4 |
 * | R | L0 | R0 | L1 | R1 | ...
 *
 *				 R
 *              / \
 *             L0 R0
 *            /  \
 *           L1  R1
 *
 */

/* Shall return true if a comes before b */
typedef bool cmpfunc(void* a, void* b);

struct prioqueue {
	u32 items;
	u32 item_size;

	struct barray* mem;

	cmpfunc* cmp;
};

static void prioqueue_heapify(struct prioqueue* pq, u32 index);

static inline void prioqueue_heapify(struct prioqueue* pq, u32 index) {
	if (pq->items <= 1)
		return;

	u32 L = PQ_LEFT(index);
	u32 R = PQ_RIGHT(index);
	u32 tmp = index;

	if (L >= pq->items || R >= pq->items)
		return;

	void* left     = barray_get(pq->mem, L);
	void* right    = barray_get(pq->mem, R);
	void* indexptr = barray_get(pq->mem, index);

	if (!pq->cmp(indexptr, left)) {
		indexptr = left;
		index = L;
	}

	if (!pq->cmp(indexptr, right)) {
		indexptr = right;
		index = R;
	}

	if (index != tmp) {
		u8 tmpbuf[pq->item_size];
		memcpy(tmpbuf, barray_get(pq->mem, tmp), pq->item_size);
		memcpy(barray_get(pq->mem, tmp), barray_get(pq->mem, index), pq->item_size);
		memcpy(barray_get(pq->mem, index), tmpbuf, pq->item_size);

		prioqueue_heapify(pq, index);
	}
}

static inline struct prioqueue* prioqueue_new(u32 items_per_bucket, u32 item_size, cmpfunc* cmp) {
	struct prioqueue* pq = (struct prioqueue*) sbmf_stack_push(sizeof(struct prioqueue));
	pq->items = 0;
	pq->item_size = item_size;
	pq->mem = barray_new(items_per_bucket, item_size);
	pq->cmp = cmp;

	return pq;
}

static inline void prioqueue_push(struct prioqueue* pq, const void* data) {
	/* Copy data to the last slot in the array */
	u32 i = pq->items;
	memcpy(barray_get(pq->mem, i), data, pq->item_size);
	pq->items++;

	void* child  = barray_get(pq->mem, i);
	void* parent = barray_get(pq->mem, PQ_PARENT(i));
	u8 tmpbuf[pq->item_size];
	while (i != 0 && pq->cmp(child, parent)) {
		/* swap child and parent */
		memcpy(tmpbuf, child, pq->item_size);
		memcpy(child, parent, pq->item_size);
		memcpy(parent, tmpbuf, pq->item_size);

		/* update indices */
		i = PQ_PARENT(i);
		child = barray_get(pq->mem, i);
		parent = barray_get(pq->mem, PQ_PARENT(i));
	}
}

static inline u8* prioqueue_top(struct prioqueue* pq) {
	return (u8*) barray_get(pq->mem, 0);
}

static inline bool prioqueue_pop(struct prioqueue* pq, void* out) {
	if (pq->items == 0)
		return false;

	memcpy(out, barray_get(pq->mem, 0), pq->item_size);

	/* swap first and last element */
	u8 tmpbuf[pq->item_size];
	memcpy(tmpbuf, barray_get(pq->mem, 0), pq->item_size);
	memcpy(barray_get(pq->mem, 0), barray_get(pq->mem, pq->items-1), pq->item_size);
	memcpy(barray_get(pq->mem, pq->items-1), tmpbuf, pq->item_size);

	pq->items--;
	prioqueue_heapify(pq, 0);

	return true;
}
