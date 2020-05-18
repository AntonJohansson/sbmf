#pragma once

#include "common.h"

// TODO: it's kinda ambiguous wheter rows,cols refer to size of entire matrix or just
// the region in data

typedef c64 mat_scalar_t;
typedef u16 mat_size_t;

typedef struct {
	bool is_row_major;
	mat_size_t rows, cols;
	mat_scalar_t* data;
} mat;

typedef struct {
	mat base;
	mat_size_t super_diags;
	mat_size_t sub_diags;
} bandmat;

typedef struct {
	mat base;
	mat_size_t bandcount;
} hermitian_bandmat;

extern void complex_hermitian_bandmat_mulv(mat_scalar_t* ans_vec, hermitian_bandmat mat, mat_scalar_t* vec);
extern void complex_bandmat_mulv(mat_scalar_t* ans_vec, bandmat mat, mat_scalar_t* vec);

#define mat_mulv(mat, vec) 																									\
	_Generic((mat),																														\
			case complex_hermitian_bandmat: hermitian_bandmat_mulv(mat, vec),			\
			case complex_bandmat: bandmat_mulv(mat, vec),													\
			default: )

extern void mat_transpose(mat* ans_mat, mat m);
