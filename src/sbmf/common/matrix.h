#pragma once

#include "common.h"

// TODO: it's kinda ambiguous wheter rows,cols refer to size of entire matrix or just
// the region in data

typedef c64 mat_scalar_t;
typedef u16 mat_size_t;

typedef struct mat {
	bool is_row_major;
	mat_size_t rows, cols;
	mat_scalar_t* data;
} mat;

typedef struct bandmat {
	mat base;
	mat_size_t super_diags;
	mat_size_t sub_diags;
} bandmat;

typedef struct hermitian_bandmat {
	mat base;
	mat_size_t bandcount;
	mat_size_t size;
} hermitian_bandmat;

mat mat_new(mat_size_t rows, mat_size_t cols);
mat mat_new_zero(mat_size_t rows, mat_size_t cols);
mat mat_duplicate(mat m);

static inline u64 mat_size(mat m) {
	return sizeof(mat_scalar_t)*m.rows*m.cols;
}

void complex_hermitian_bandmat_mulv(mat_scalar_t* ans_vec, hermitian_bandmat mat, mat_scalar_t* vec);
void complex_bandmat_mulv(mat_scalar_t* ans_vec, bandmat mat, mat_scalar_t* vec);

#define mat_mulv(mat, vec) 																									\
	_Generic((mat),																														\
			case complex_hermitian_bandmat: hermitian_bandmat_mulv(mat, vec),			\
			case complex_bandmat: bandmat_mulv(mat, vec),													\
			default: )

void mat_transpose_raw(mat_scalar_t* ans, mat_scalar_t* in, mat_size_t rows, mat_size_t cols);
void mat_transpose(mat* ans_mat, mat m);

hermitian_bandmat construct_finite_diff_mat(u32 samples_per_dimension, u32 dimensions, f64* deltas);
