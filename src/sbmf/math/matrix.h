#pragma once

#include <sbmf/sbmf.h>
#include <sbmf/types.h>
#include <string.h>

/* TODO: it's kinda ambiguous wheter rows,cols refer to size of entire matrix or just
 * the region in data
 */

struct complex_hermitian_bandmat {
	c64* data;
	u32 bandcount; /* rows */
	u32 size; /* cols */
};

static inline struct complex_hermitian_bandmat complex_hermitian_bandmat_new(u32 bandcount, u32 size) {
	return (struct complex_hermitian_bandmat) {
		.data = (c64*)sbmf_stack_push(sizeof(c64)*size*bandcount),
		.bandcount = bandcount,
		.size = size,
	};
}

static inline struct complex_hermitian_bandmat complex_hermitian_bandmat_new_zero(u32 bandcount, u32 size) {
	struct complex_hermitian_bandmat bm = complex_hermitian_bandmat_new(bandcount, size);
	memset(bm.data, 0, sizeof(c64)*size*bandcount);
	return bm;
}

/*
 *
 *				(0,0)	(0,1)
 *						(1,1)	(1,2)
 *								(2,2)	(2,3)
 *										(3,3)
 *
 *				(0,0)	(0,1)	(0,2)	(0,3)
 *						(1,1)	(1,2)	(1,3)
 *								(2,2)	(2,3)
 *										(3,3)
 *
 * 				for r = 0 -> 4
 * 					for c = r -> 4
 *
 *
 */

#define U32MIN(a,b) \
	((a < b) ? a : b)

#define COMPLEX_HERMITIAN_BANDMAT_FOREACH(bm, r,c) 						\
	for (u32 r = 0; r < bm.size; ++r)									\
		for (u32 c = r; c < U32MIN(bm.size, r+bm.bandcount); ++c)

/* How do we compute the position in band storage?
 *
 * 			bandcount == size == 4
 *
 *		(x,x) (x,x) (x,x) (0,3)		0	1	2	3
 *		(x,x) (x,x) (0,2) (1,3)		4	5	6	7
 *		(x,x) (0,1) (1,2) (2,3)		8	9	10	11
 *		(0,0) (1,1) (2,2) (3,3)		12	13	14	15
 *
 * 			bandcount = 3, size = 4
 *
 *		(x,x) (x,x) (0,2) (1,3)		0	1	2	3
 *		(x,x) (0,1) (1,2) (2,3)		4	5	6 	7
 *		(0,0) (1,1) (2,2) (3,3)		8 	9 	10	11
 *
 *		r,c -> (n-1)*(n - (c-r)) + r - (n-b)*n
 *		       n(n-1) - (c-r)(n-1) + r - n(n-b)
 *		       n(n-1-(n-b)) - (c-r)(n-1) + r
 *		       n(b-1) - (c-r)(n-1) + r
 *		       n(b-1) - c(n-1) + r(n-1) + r
 *		       n(b-1) - c(n-1) + rn
 *		       n(b-1) - cn + c + rn
 *		       n(b-1+(r-c)) + c
 *
 *		(x,x) (0,1) (1,2) (2,3)		0	1	2 	3
 *		(0,0) (1,1) (2,2) (3,3)		4 	5 	6 	7
 *
 * 				bandcount, column coord
 *
 *		(x,x) (1,1) (1,2) (1,3)		0	1	2 	3
 *		(0,0) (0,1) (0,2) (0,3)		4 	5 	6 	7
 *
 * 				bandcount, column coord
 *
 *		(x,x) (1,1) (1,2) (1,3)		0	1	2 	3
 *		(0,0) (0,1) (0,2) (0,3)		4 	5 	6 	7
 *
 *		(x,x) (0,0)		0	4            0 1
 *		(0,1) (1,1)		1	5            2 3
 *		(1,2) (2,2)		2	6     ->     4 5
 *		(2,3) (3,3)		3	7            6 7
 *
 *		r,c -> n(b-1+(r-c)) + c
 *
 *		(0,1) -> 4*(2 - 1 + (0-1)) + 1 = 1
 *		(1,0) -> 4*(2 - 1 + (1-0)) + 0 = 1
 *
 *		(1,0) -> 2*(4 - 1 + (1-0)) + 0 = 5
 *
 *
 *
 */

static inline u32 complex_hermitian_bandmat_index(struct complex_hermitian_bandmat bm, u32 row, u32 col) {
	return bm.size * (bm.bandcount - 1 + (row - col)) + col;
}

void complex_hermitian_bandmat_mulv(c64* ans_vec, struct complex_hermitian_bandmat bm, c64* vec);

struct complex_hermitian_bandmat construct_finite_diff_mat(u32 samples_per_dimension, u32 dimensions, f64* deltas);

static inline bool complex_hermitian_bandmat_is_valid(struct complex_hermitian_bandmat bm) {
	for (u32 r = 0; r < bm.size; ++r) {
		for (u32 c = r; c < bm.bandcount; ++c) {
			u32 i = complex_hermitian_bandmat_index(bm, r,c);
			if (!f64_is_valid(creal(bm.data[i])) || !f64_is_valid(cimag(bm.data[i])))
				return false;
		}
	}
	return true;
}
