#include <sbmf/math/matrix.h>
#include <sbmf/sbmf.h>

#include <cblas.h>

#include <assert.h>
#include <string.h> // memcpy, memset
#include <stdio.h> // sprintf

void hermitian_bandmat_mulv(f64* ans_vec, struct hermitian_bandmat bm, f64* vec) {
	static const f64 one = 1, zero = 0;

	const u32 num_super_diags = bm.bandcount-1;

	f64 bmtrans[bm.size*bm.bandcount];
	for (u32 r = 0; r < bm.bandcount; ++r) {
		for (u32 c = 0; c < bm.size; ++c) {
			bmtrans[c*bm.bandcount + r] = bm.data[r*bm.size + c];
		}
	}

	cblas_dsbmv(CblasColMajor, CblasUpper,
			bm.size, num_super_diags,
			one, bmtrans, bm.bandcount, vec, 1, zero, ans_vec, 1);
}

void complex_hermitian_bandmat_mulv(c64* ans_vec, struct complex_hermitian_bandmat bm, c64* vec) {
	static const c64 one = 1, zero = 0;

	const u32 num_super_diags = bm.bandcount-1;

	c64 bmtrans[bm.size*bm.bandcount];
	for (u32 r = 0; r < bm.bandcount; ++r) {
		for (u32 c = 0; c < bm.size; ++c) {
			bmtrans[c*bm.bandcount + r] = bm.data[r*bm.size + c];
		}
	}

	cblas_zhbmv(CblasColMajor, CblasUpper,
			bm.size, num_super_diags,
			&one, bmtrans, bm.bandcount, vec, 1, &zero, ans_vec, 1);

}

struct complex_hermitian_bandmat construct_finite_diff_mat(u32 samples_per_dimension, u32 dimensions, f64* deltas) {
	u32 size = pow(samples_per_dimension, dimensions);
	u32 bands = pow(samples_per_dimension, dimensions-1);
	struct complex_hermitian_bandmat m = complex_hermitian_bandmat_new_zero(bands+1, size);

	// Setup main diagonal
	f64 main_diag_value = pow(2, dimensions);
	for (u32 i = (m.bandcount-1)*size; i < (m.bandcount)*size; ++i) {
		m.data[i] = -main_diag_value/(deltas[0]*deltas[0]);
	}

	// Setup off-diagonal elements
	for (u32 i = 1; i <= dimensions; ++i) {
		u32 bandindex = pow(samples_per_dimension, i-1);
		for (u32 j = 0; j < size; ++j) {
			m.data[j + (m.bandcount-1 - bandindex)*size] = 1.0/(deltas[i-1]*deltas[i-1]);
		}
	}

	return m;
}

void complex_hermitian_bandmat_print(struct complex_hermitian_bandmat bm, const char label[]) {
	sbmf_log_info("[%ux%u] %s", bm.size, bm.size, label);
	COMPLEX_HERMITIAN_BANDMAT_FOREACH(bm, r,c) {
		u32 index = complex_hermitian_bandmat_index(bm, r,c);
		printf("%.2e+%.2ei\t", CCOMP(bm.data[index]));
		if (c == bm.size-1) {
			printf("\n");
		}
	}
	printf("\n");
}
