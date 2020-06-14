#include "matrix.h"
#include <string.h> // memcpy, memset
#include <cblas.h>
#include <sbmf/memory/stack_allocator.h>
#include <assert.h>

mat mat_new(mat_size_t rows, mat_size_t cols) {
	return (mat){
		.is_row_major = true,
		.rows = rows,
		.cols = cols,
		.data = (mat_scalar_t*)sa_push(_sbmf.main_stack, sizeof(mat_scalar_t)*rows*cols)
	};
}

mat mat_new_zero(mat_size_t rows, mat_size_t cols) {
	mat m = mat_new(rows, cols);
	memset(m.data, 0, mat_size(m));
	return m;
}

mat mat_duplicate(mat m) {
	mat n = mat_new(m.rows, m.cols);
	n.is_row_major = m.is_row_major;
	return n;
}

void complex_hermitian_bandmat_mulv(mat_scalar_t* ans_vec, hermitian_bandmat bmat, mat_scalar_t* vec) {
	static const mat_scalar_t one = 1, zero = 0;
	const mat_size_t num_super_diags = bmat.bandcount-1;
	// TODO: Passing data with is_row_major doesnt produce correct result ://
	// manually transpose?
	if (bmat.base.is_row_major) {
		mat tmp = mat_duplicate(bmat.base);
		mat_transpose(&tmp, bmat.base);
		cblas_zhbmv(CblasColMajor, CblasUpper,
				bmat.size, num_super_diags,
				&one, tmp.data, bmat.bandcount, vec, 1, &zero, ans_vec, 1);
	} else {
		cblas_zhbmv(CblasColMajor, CblasUpper,
				bmat.size, num_super_diags,
				&one, bmat.base.data, bmat.bandcount, vec, 1, &zero, ans_vec, 1);
	}

}

void complex_bandmat_mulv(mat_scalar_t* ans_vec, bandmat mat, mat_scalar_t* vec) {
	static const mat_scalar_t one = 1, zero = 0;
	cblas_zgbmv((mat.base.is_row_major) ? CblasRowMajor : CblasColMajor, CblasNoTrans,
			mat.base.rows, mat.base.cols, mat.sub_diags, mat.super_diags,
			&one, mat.base.data, mat.base.rows, vec, 1, &zero, ans_vec, 1);
}

void mat_transpose_raw(mat_scalar_t* ans, mat_scalar_t* in, mat_size_t rows, mat_size_t cols) {
	// TODO: this can be freed
	mat_scalar_t* temp = (mat_scalar_t*)sa_push(_sbmf.main_stack, sizeof(mat_scalar_t)*rows*cols);

	for (mat_size_t r = 0; r < rows; ++r) {
		for (mat_size_t c = 0; c < cols; ++c) {
			temp[r + c*(rows)] = in[c + r*(cols)];
		}
	}

	// rows-1 + (cols-1)*rows
	// rows-1 + cols*rows - rows
	// cols*rows-1

	memcpy(ans, temp, sizeof(mat_scalar_t)*rows*cols);
}

void mat_transpose(mat* ans_mat, mat m) {
	assert(m.rows*m.cols == ans_mat->rows*ans_mat->cols);
	mat_transpose_raw(ans_mat->data, m.data, m.rows, m.cols);

	ans_mat->is_row_major = !m.is_row_major;
	mat_size_t swap_temp = m.rows;
	ans_mat->rows = m.cols;
	ans_mat->cols = swap_temp;
}

hermitian_bandmat construct_finite_diff_mat(u32 samples_per_dimension, u32 dimensions, f64* deltas) {
	i32 size = pow(samples_per_dimension, dimensions);
	i32 bands = pow(samples_per_dimension, dimensions-1);
	hermitian_bandmat m = {
		.base = mat_new_zero(bands+1,size),
		.bandcount = bands+1,
		.size = size,
	};

	// Setup main diagonal
	i32 main_diag_value = pow(2, dimensions);
	for (i32 i = (m.bandcount-1)*size; i < (m.bandcount)*size; ++i) {
		m.base.data[i] = -main_diag_value/(deltas[0]*deltas[0]);
	}

	// Setup off-diagonal elements
	for (i32 i = 1; i <= dimensions; ++i) {
		i32 bandindex = pow(samples_per_dimension, i-1);
		for (i32 j = 0; j < size; ++j) {
			m.base.data[j + (m.bandcount-1 - bandindex)*size] = 1/(deltas[i-1]*deltas[i-1]);
		}
	}

	return m;
}
