f64 gaussian(f64 x, f64 mu, f64 sigma) {
	return 1.0/(sigma*sqrt(2.0*M_PI)) * exp(-(x-mu)*(x-mu)/(2*sigma*sigma));
}

void f64_normalize(f64* out, f64* data, u32 size) {
	f64 sum = 0.0;

#pragma omp parallel for shared(data) reduction(+: sum)
	for (u32 i = 0; i < size; ++i) {
		f64 absval = fabs(data[i]);
		sum += absval*absval;
	}

	f64 scaling = 1.0/sqrt(sum);
#pragma omp parallel for
	for (u32 i = 0; i < size; ++i) {
		out[i] = data[i] * scaling;
	}
}

static inline f128 factorial_128(const u32 n) {
	f128 prod = 1.0;
	f128 current_value = (f128) n;
	while (current_value > 0.0) {
		prod *= current_value;
		current_value -= 1.0;
	}
	return prod;
}

static inline u64 n_choose_k(const u32 n, const u32 k) {
	f128 n_fact = factorial_128(n);
	f128 k_fact = factorial_128(k);
	f128 n_minus_k_fact = factorial_128(n-k);
	f128 ans = (n_fact/n_minus_k_fact) * 1.0/k_fact;
	return lroundl(ans);
}
