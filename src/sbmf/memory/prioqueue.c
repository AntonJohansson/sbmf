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
	u32 max_items;

	//struct barray* mem;
	void* mem;

	cmpfunc* cmp;
};

static void* array_get(struct prioqueue* pq, u32 index) {
	return (u8*)pq->mem + index*pq->item_size;
}

static void prioqueue_heapify(struct prioqueue* pq, u32 index);

static inline void prioqueue_heapify(struct prioqueue* pq, u32 index) {
	if (pq->items <= 1)
		return;

	u32 L = PQ_LEFT(index);
	u32 R = PQ_RIGHT(index);
	u32 tmp = index;

	void* left     = array_get(pq, L);
	void* right    = array_get(pq, R);
	void* indexptr = array_get(pq, index);

	if (L < pq->items) {
		if (!pq->cmp(indexptr, left)) {
			indexptr = left;
			index = L;
		}
	}

	if (R < pq->items) {
		if (!pq->cmp(indexptr, right)) {
			indexptr = right;
			index = R;
		}
	}

	if (index != tmp) {
		u8 tmpbuf[pq->item_size];
		memcpy(tmpbuf, array_get(pq, tmp), pq->item_size);
		memcpy(array_get(pq, tmp), array_get(pq, index), pq->item_size);
		memcpy(array_get(pq, index), tmpbuf, pq->item_size);

		prioqueue_heapify(pq, index);
	}
}

static inline struct prioqueue* prioqueue_new(void* mem, u32 items_per_bucket, u32 item_size, cmpfunc* cmp) {
	struct prioqueue* pq = (struct prioqueue*) sbmf_stack_push(sizeof(struct prioqueue));
	pq->items = 0;
	pq->item_size = item_size;
	pq->max_items = items_per_bucket;
	//pq->mem = barray_new(items_per_bucket, item_size);
	pq->mem = mem;
	pq->cmp = cmp;

	return pq;
}

static inline void prioqueue_push(struct prioqueue* pq, const void* data) {
	assert(pq->items < pq->max_items);

	/* Copy data to the last slot in the array */
	u32 i = pq->items;
	memcpy(array_get(pq, i), data, pq->item_size);
	pq->items++;

	void* child  = array_get(pq, i);
	void* parent = array_get(pq, PQ_PARENT(i));
	u8 tmpbuf[pq->item_size];
	/* While the child comes before the parent */
	while (i != 0 && pq->cmp(child, parent)) {
		/* swap child and parent */
		memcpy(tmpbuf, child, pq->item_size);
		memcpy(child, parent, pq->item_size);
		memcpy(parent, tmpbuf, pq->item_size);

		/* update indices */
		i = PQ_PARENT(i);
		child = array_get(pq, i);
		parent = array_get(pq, PQ_PARENT(i));
	}
}

static inline u8* prioqueue_top(struct prioqueue* pq) {
	return (u8*) array_get(pq, 0);
}

static inline bool prioqueue_pop(struct prioqueue* pq, void* out) {
	if (pq->items == 0)
		return false;

	if (out)
		memcpy(out, array_get(pq, 0), pq->item_size);

	/* swap first and last element */
	u8 tmpbuf[pq->item_size];
	memcpy(tmpbuf, array_get(pq, 0), pq->item_size);
	memcpy(array_get(pq, 0), array_get(pq, pq->items-1), pq->item_size);
	memcpy(array_get(pq, pq->items-1), tmpbuf, pq->item_size);

	pq->items--;
	prioqueue_heapify(pq, 0);

	return true;
}
