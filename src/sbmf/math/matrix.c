struct hermitian_bandmat hermitian_bandmat_new(u32 bandcount, u32 size) {
	return (struct hermitian_bandmat) {
		.data = (f64*)sbmf_stack_push(sizeof(f64)*size*bandcount),
		.bandcount = bandcount,
		.size = size,
	};
}

struct hermitian_bandmat hermitian_bandmat_new_zero(u32 bandcount, u32 size) {
	struct hermitian_bandmat bm = hermitian_bandmat_new(bandcount, size);
	memset(bm.data, 0, sizeof(f64)*size*bandcount);
	return bm;
}

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

static bool hermitian_bandmat_is_valid(struct hermitian_bandmat bm) {
	//f64 smallest_abs =  INFINITY;
	//f64 largest_abs  = -INFINITY;
	HERMITIAN_BANDMAT_FOREACH(bm, r,c) {
		u32 i = hermitian_bandmat_index(bm, r,c);
		if (isnan(bm.data[i]) || isinf(bm.data[i]))
			return false;
	}

	return true;
}
