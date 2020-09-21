#pragma once

#include <sbmf/sbmf.h>
#include <string.h> /* memcpy, memove */

#define PQ_PARENT(x) ((x-1)/2)
#define PQ_LEFT(x)	 (2*x+1)
#define PQ_RIGHT(x)  (2*x+2)
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

struct pqheap {
	u32 items;
	u32 max_items;
	u32 item_size;

	u8* mem;

	cmpfunc* cmp;
};

void pqheap_heapify(struct pqheap* pq, u32 index);

inline void pqheap_heapify(struct pqheap* pq, u32 index) {
	if (pq->items <= 1)
		return;

	u32 L = PQ_LEFT(index);
	u32 R = PQ_RIGHT(index);
	u32 tmp = index;

	if (L >= pq->items || R >= pq->items)
		return;

	void* left     = &pq->mem[L     * pq->item_size];
	void* right    = &pq->mem[R     * pq->item_size];
	void* indexptr = &pq->mem[index * pq->item_size];

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
		memcpy(tmpbuf, &pq->mem[tmp * pq->item_size], pq->item_size);
		memcpy(&pq->mem[tmp * pq->item_size], &pq->mem[index * pq->item_size], pq->item_size);
		memcpy(&pq->mem[index * pq->item_size], tmpbuf, pq->item_size);

		pqheap_heapify(pq, index);
	}
}

inline struct pqheap* pqheap_new(u32 items, u32 item_size, cmpfunc* cmp) {
	u8* mem = sbmf_stack_push(items*item_size + sizeof(struct pqheap));

	struct pqheap* pq = (struct pqheap*) mem;
	pq->max_items = items;
	pq->item_size = item_size;
	pq->mem = (u8*) (pq + 1);
	pq->cmp = cmp;

	return pq;
}

inline void pqheap_push(struct pqheap* pq, const void* data) {
	/* Copy data to the last slot in the array */
	u32 i = pq->items;
	memcpy(&pq->mem[i * pq->item_size], data, pq->item_size);
	pq->items++;

	void* child  = &pq->mem[i            * pq->item_size];
	void* parent = &pq->mem[PQ_PARENT(i) * pq->item_size];
	u8 tmpbuf[pq->item_size];
	while (i != 0 && pq->cmp(child, parent)) {
		/* swap child and parent */
		memcpy(tmpbuf, child, pq->item_size);
		memcpy(child, parent, pq->item_size);
		memcpy(parent, tmpbuf, pq->item_size);

		/* update indices */
		i = PQ_PARENT(i);
		child = &pq->mem[i             * pq->item_size];
		parent = &pq->mem[PQ_PARENT(i) * pq->item_size];
	}
}

inline u8* pqheap_top(struct pqheap* pq) {
	return &pq->mem[0];
}

inline bool pqheap_pop(struct pqheap* pq) {
	if (pq->items == 0)
		return false;

	/* swap first and last element */
	u8 tmpbuf[pq->item_size];
	memcpy(tmpbuf, &pq->mem[0], pq->item_size);
	memcpy(&pq->mem[0], &pq->mem[(pq->items-1)*pq->item_size], pq->item_size);
	memcpy(&pq->mem[(pq->items-1)*pq->item_size], tmpbuf, pq->item_size);

	pq->items--;
	pqheap_heapify(pq, 0);

	return true;
}
