f64 hermite_integral_4(u32 i, u32 j, u32 k, u32 l) {
	if ((i+j+k+l) % 2 != 0)
		return 0.0;

	u32 m_max = (i < j) ? i : j;

	f128 sum = 0.0;
	for (u32 m = 0; m <= m_max; ++m) {

		f128 df1 = double_factorial_128( (i+j-2*m) + k - l - 1);
		f128 df2 = double_factorial_128( (i+j-2*m) - k + l - 1);
		f128 df3 = double_factorial_128(-(i+j-2*m) + k + l - 1);

		f128 f1 = factorial_128(m);
		f128 f2 = factorial_128(i-m);
		f128 f3 = factorial_128(j-m);

		sum += (pow(2,m)/f1)*(1/f2)*(1/f3)*(df1*df2*df3);
	}

	const f64 sqrt_omega_over_2pi = sqrt(OMEGA/(2.0*M_PI));
	f128 f1 = factorial_128(i);
	f128 f2 = factorial_128(j);
	f128 f3 = factorial_128(k);
	f128 f4 = factorial_128(l);

	return sqrt_omega_over_2pi*sqrt(((f1/f3)*(f2/f4)))/pow(2,(i+j+k+l)/2)*sum;
}
