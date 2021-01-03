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

inline u32 diis_get(struct diis_data* d, const u32 log_entry, const u32 index_in_entry) {
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
	/* Compute deltas */
	for (u32 i = 0; i < d->log_length-1; ++i) {
		u32 l0 = (d->log_index + i + 0) % d->log_length;
		u32 l1 = (d->log_index + i + 1) % d->log_length;
		//printf("%u,%u\n", l0,l1);
		for (u32 j = 0; j < d->size; ++j) {
			*(d->delta_log + diis_get(d, i,j)) =
				*(d->data_log + diis_get(d, l1,j)) - *(d->data_log + diis_get(d, l0,j));
		}
	}


	/* Construct matrix */
	f64 mat[d->log_length * d->log_length];
	{
		for (u32 r = 0; r < d->log_length-1; ++r) {
			for (u32 c = r; c < d->log_length-1; ++c) {
				/* Compute overlap */
				f64 overlap = 0.0;
				for (u32 i = 0; i < d->size; ++i) {
					f64 a = *(d->delta_log + diis_get(d, r,i));
					f64 b = *(d->delta_log + diis_get(d, c,i));
					overlap += a*b;
				}

				mat[r*d->log_length + c] = overlap;
				if (r != c) {
					mat[c*d->log_length + r] = overlap;
				}
			}
		}

		/* Set bottom and right edges to -1 */
		for (u32 i = 0; i < d->log_length-1; ++i) {
			u32 r = d->log_length-1;
			mat[r*d->log_length + i] = -1;
			mat[i*d->log_length + r] = -1;
		}

		/* Set bottom right to 0 */
		mat[(d->log_length-1)*d->log_length + (d->log_length-1)] = 0;
	}

	for (u32 r = 0; r < d->log_length; ++r) {
		for (u32 c = 0; c < d->log_length; ++c) {
			printf("%lf ", mat[r*d->log_length + c]);
		}
		printf("\n");
	}

	/* Solve mat*x = y problem */
	f64 y[d->log_length];
	memset(y, 0, d->log_length*sizeof(f64));
	y[d->log_length-1] = -1;

	i32 pivot[d->log_length];

	/* LU factorize mat */
	{
		LAPACKE_dgetrf(LAPACK_COL_MAJOR, d->log_length, d->log_length,
				mat, d->log_length, pivot);
	}

	/* Solve problem with LU factorized mat */
	{
		LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', d->log_length,
				1, mat, d->log_length,
				pivot,
				y, d->log_length);
	}

	memset(new_x, 0, d->size*sizeof(f64));
	*new_E = 0;
	/* optimal coeffs. are now stored in y */
	for (u32 i = 0; i < d->log_length-1; ++i) {
		u32 l = (d->log_index + i) % d->log_length;
		if (isnan(y[i]) || isinf(y[i]))
			return false;

		*new_E += y[i]*d->energy[i];

		for (u32 j = 0; j < d->size; ++j) {
			new_x[j] += y[i] * (*(d->data_log + diis_get(d, l,j)));
		}
	}

	return true;
}
