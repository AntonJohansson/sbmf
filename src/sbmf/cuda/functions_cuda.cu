__host__
f64 gaussian(f64 x, f64 mu, f64 sigma);

__host__
void f64_normalize(f64* out, f64* data, u32 size);

__host__
static inline f128 factorial_128_cuda(const u32 n) {
	f128 prod = 1.0;
	f128 current_value = (f128) n;
	while (current_value > 0.0) {
		prod *= current_value;
		current_value -= 1.0;
	}
	return prod;
}

__host__
static inline f128 double_factorial_positive_128_cuda(const u32 n) {
	f128 prod = 1.0;
	f128 current_value = (f128) n;
	while (current_value > 0.0) {
		prod *= current_value;
		current_value -= 2.0;
	}
	return prod;
}

__host__
static inline f128 double_factorial_128_cuda(const i32 n) {
	if (n < 0) {
		assert(n % 2 != 0);
		f128 n_double_fact = double_factorial_positive_128_cuda((u32)(-n));
		f128 sign = (((n-1)/2) % 2 == 0) ? 1.0 : -1.0;
		return sign*n/n_double_fact;
	}

	return double_factorial_positive_128_cuda((u32)n);
}

__host__
static inline u64 n_choose_k_cuda(const u32 n, const u32 k) {
	f128 n_fact = factorial_128_cuda(n);
	f128 k_fact = factorial_128_cuda(k);
	f128 n_minus_k_fact = factorial_128_cuda(n-k);
	f128 ans = (n_fact/n_minus_k_fact) * 1.0/k_fact;
	return lroundl(ans);
}
