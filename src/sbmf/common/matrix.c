#include "matrix.h"
#include <string.h> // memcpy
#include <cblas.h>

void complex_hermitian_bandmat_mulv(mat_scalar_t* ans_vec, hermitian_bandmat mat, mat_scalar_t* vec) {
	static const mat_scalar_t one = 1, zero = 0;
	const mat_size_t num_super_diags = mat.bandcount-1;
	cblas_zhbmv((mat.base.is_row_major) ? CblasRowMajor : CblasColMajor, CblasUpper, 
			mat.base.rows*mat.base.cols, num_super_diags,
			&one, mat.base.data, mat.base.rows, vec, 1, &zero, ans_vec, 1);
}

void complex_bandmat_mulv(mat_scalar_t* ans_vec, bandmat mat, mat_scalar_t* vec) {
	static const mat_scalar_t one = 1, zero = 0;
	cblas_zgbmv((mat.base.is_row_major) ? CblasRowMajor : CblasColMajor, CblasNoTrans, 
			mat.base.rows, mat.base.cols, mat.sub_diags, mat.super_diags, 
			&one, mat.base.data, mat.base.rows, vec, 1, &zero, ans_vec, 1);
}

void mat_transpose(mat* ans_mat, mat m) {
	const mat_size_t size = m.rows*m.cols;
	mat_scalar_t temp[size];

	for (mat_size_t c = 0; c < m.cols; ++c) {
		for (mat_size_t r = 0; r < m.rows; ++r) {
			temp[r + c*m.rows] = m.data[c + r*m.cols];
		}
	}

	mat_size_t swap_temp = ans_mat->rows;
	ans_mat->rows = ans_mat->cols;
	ans_mat->cols = swap_temp;
	ans_mat->is_row_major = !m.is_row_major;

	memcpy(ans_mat->data, temp, sizeof(mat_scalar_t)*size);
}
