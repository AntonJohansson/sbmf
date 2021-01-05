static bool symmetric_bandmat_is_valid(struct symmetric_bandmat bm) {
	//f64 smallest_abs =  INFINITY;
	//f64 largest_abs  = -INFINITY;
	SYMMETRIC_BANDMAT_FOREACH(bm, r,c) {
		u32 i = symmetric_bandmat_index(bm, r,c);
		if (isnan(bm.data[i]) || isinf(bm.data[i]))
			return false;
	}

	return true;
}

struct symmetric_bandmat symmetric_bandmat_new(u32 bandcount, u32 size) {
	return (struct symmetric_bandmat) {
		.data 		= sbmf_stack_push(sizeof(f64)*size*bandcount),
		.bandcount 	= bandcount,
		.size 		= size,
	};
}

struct symmetric_bandmat symmetric_bandmat_new_zero(u32 bandcount, u32 size) {
	struct symmetric_bandmat bm = symmetric_bandmat_new(bandcount, size);
	memset(bm.data, 0, sizeof(f64)*size*bandcount);
	return bm;
}

void symmetric_bandmat_mulv(f64* ans_vec, struct symmetric_bandmat bm, f64* vec) {
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
