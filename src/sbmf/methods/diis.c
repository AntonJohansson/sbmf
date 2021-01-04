struct diis_data {
	f64* data_log;
	f64* delta_log;
	f64* energy;
	u32 log_length;
	u32 size;

	u32 log_index;
	u32 log_items;
};

struct diis_data* diis_data_new(const u32 size, const u32 log_length) {
	void* mem = sbmf_stack_push(
			sizeof(struct diis_data)
			+ size*log_length*sizeof(f64) /* data log */
			+ size*(log_length-1)*sizeof(f64) /* delta log */
			+ log_length*sizeof(f64) /* energy */
			);

	struct diis_data* d = mem;
	d->data_log  = (f64*)((u8*)mem + sizeof(struct diis_data));
	d->delta_log = (f64*)((u8*)mem + sizeof(struct diis_data) + size*log_length*sizeof(f64));
	d->energy 	 = (f64*)((u8*)mem + sizeof(struct diis_data) + size*log_length*sizeof(f64) + size*(log_length-1)*sizeof(f64));


	d->log_length = log_length;
	d->size = size;
	d->log_index = 0;
	d->log_items = 0;

	return d;
}

static inline u32 diis_get(struct diis_data* d, const u32 log_entry, const u32 index_in_entry) {
	return log_entry * d->size + index_in_entry;
}

void diis_push(struct diis_data* d, f64* v, f64 E) {
	memcpy(d->data_log + diis_get(d, d->log_index, 0), v, d->size * sizeof(f64));
	d->energy[d->log_index] = E;
	d->log_index = (d->log_index + 1) % d->log_length;

	if (d->log_items < d->log_length)
		d->log_items++;
}

bool diis_compute(struct diis_data* d, f64* new_x, f64* new_E) {
	if (d->log_items < 3)
		return false;

	/* Compute deltas */
	for (u32 i = 0; i < d->log_items-1; ++i) {
		u32 l0 = (d->log_index + i + 0) % d->log_items;
		u32 l1 = (d->log_index + i + 1) % d->log_items;
		//printf("%u,%u\n", l0,l1);
		for (u32 j = 0; j < d->size; ++j) {
			*(d->delta_log + diis_get(d, i,j)) =
				*(d->data_log + diis_get(d, l1,j)) - *(d->data_log + diis_get(d, l0,j));
		}
	}

	/* Construct matrix */
	f64 mat[d->log_items * d->log_items];
	{
		for (u32 r = 0; r < d->log_items-1; ++r) {
			for (u32 c = r; c < d->log_items-1; ++c) {
				/* Compute overlap */
				f64 overlap = 0.0;
				for (u32 i = 0; i < d->size; ++i) {
					f64 a = *(d->delta_log + diis_get(d, r,i));
					f64 b = *(d->delta_log + diis_get(d, c,i));
					overlap += a*b;
				}

				mat[r*d->log_items + c] = overlap;
				if (r != c) {
					mat[c*d->log_items + r] = overlap;
				}
			}
		}

		/* Set bottom and right edges to -1 */
		for (u32 i = 0; i < d->log_items-1; ++i) {
			u32 r = d->log_items-1;
			mat[r*d->log_items + i] = -1;
			mat[i*d->log_items + r] = -1;
		}

		/* Set bottom right to 0 */
		mat[(d->log_items-1)*d->log_items + (d->log_items-1)] = 0;
	}

	//for (u32 r = 0; r < d->log_items; ++r) {
	//	for (u32 c = 0; c < d->log_items; ++c) {
	//		printf("%lf ", mat[r*d->log_items + c]);
	//	}
	//	printf("\n");
	//}

	/* Solve mat*x = y problem */
	f64 y[d->log_items];
	memset(y, 0, d->log_items*sizeof(f64));
	y[d->log_items-1] = -1;

	i32 pivot[d->log_items];

	/* LU factorize mat */
	{
		LAPACKE_dgetrf(LAPACK_COL_MAJOR, d->log_items, d->log_items,
				mat, d->log_items, pivot);
	}

	/* Solve problem with LU factorized mat */
	{
		LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', d->log_items,
				1, mat, d->log_items,
				pivot,
				y, d->log_items);
	}

	f64 delta[d->size];
	memset(delta, 0, d->size*sizeof(f64));
	for (u32 i = 0; i < d->log_items-1; ++i) {
		u32 l = (d->log_index + i) % d->log_items;
		for (u32 j = 0; j < d->size; ++j) {
			delta[j] += y[i] * (*(d->delta_log + diis_get(d, l,j)));
		}
	}

	/*
	 * The goal of this whole things is to minimize the following
	 * in a least squared sense to approximate the null vector
	 */
	f64 sum = 0;
	//for (u32 j = 0; j < d->size; ++j) {
	//	sum += delta[j]*delta[j];
	//}
	for (u32 r = 0; r < d->log_items-1; ++r) {
		for (u32 c = 0; c < d->log_items-1; ++c) {
			/* Compute overlap */
			f64 overlap = 0.0;
			for (u32 i = 0; i < d->size; ++i) {
				f64 a = *(d->delta_log + diis_get(d, r,i));
				f64 b = *(d->delta_log + diis_get(d, c,i));
				overlap += a*b;
			}

			sum += y[r]*y[c]*overlap;
		}
	}
	sum = sqrt(sum);
	printf("\t\tsum: %e\n", sum);

	/*
	 * When close to convergence, some equation in mat
	 * may be linearly dependant, thus some y[i] may be
	 * inf or nan. When this happens, empty the subspace.
	 */
	if (isnan(sum) || isinf(sum)) {
		d->log_items = 0;
		return false;
	}

	//for (u32 j = 0; j < d->size; ++j) {
	//	printf("%lf ", delta[j]);
	//}
	//printf("\n");

	//for (u32 i = 0; i < d->log_items; ++i) {
	//	printf("[%u]: %lf\n", i, d->energy[i]);
	//}

	memset(new_x, 0, d->size*sizeof(f64));
	*new_E = 0;
	/* optimal coeffs. are now stored in y */
	for (u32 i = 0; i < d->log_items-1; ++i) {
		u32 l = (d->log_index + i) % d->log_items;
		if (isnan(y[i]) || isinf(y[i]))
			return false;

		*new_E += y[i]*d->energy[i];

		for (u32 j = 0; j < d->size; ++j) {
			new_x[j] += y[i] * (*(d->data_log + diis_get(d, l,j)));
		}
	}

	return true;
}
