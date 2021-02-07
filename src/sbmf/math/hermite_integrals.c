#define SQRT_PI_OVER_2 (1.2533141373155002512078826424055226265034933703049691583149617881)

static inline f128 hermite_integral_3(u32 i, u32 j, u32 k) {
	if ((i+j+k) % 2 != 0)
		return 0.0;

	f128 f1 = double_factorial_128( i + j - k - 1);
	f128 f2 = double_factorial_128( i - j + k - 1);
	f128 f3 = double_factorial_128(-i + j + k - 1);
	return f1*f2*f3;
}

static inline f64 hermite_integral_4(u32 i, u32 j, u32 k, u32 l) {
	u32 m_max = (i < j) ? i : j;

	f128 sum = 0.0;
	for (u32 m = 0; m <= m_max; ++m) {
		f128 integral = hermite_integral_3(i+j-2*m, k, l);
		if (integral == 0)
			continue;

		f128 b1 = n_choose_k(i, m);
		f128 b2 = n_choose_k(j, m);
		f128 m_fact = factorial_128(m);

		sum += (ho_K(i)*b1)*(ho_K(j)*b2)*(ho_K(k)*m_fact)*(ho_K(l)*pow(2,m))*integral;
	}

	return SQRT_PI_OVER_2*sum;
}
