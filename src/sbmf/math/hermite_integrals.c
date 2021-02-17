#include <fenv.h>

f64 hermite_integral_4(u32 i, u32 j, u32 k, u32 l) {
	if ((i+j+k+l) % 2 != 0)
		return 0.0;

	feclearexcept(FE_ALL_EXCEPT);

	u32 m_max = (i < j) ? i : j;

	f128 ni = 1.0/sqrt(factorial_128(i)*pow(2,i));
	f128 nj = 1.0/sqrt(factorial_128(j)*pow(2,j));
	f128 nk = 1.0/sqrt(factorial_128(k)*pow(2,k));
	f128 nl = 1.0/sqrt(factorial_128(l)*pow(2,l));

	f128 sum = 0.0;
	for (u32 m = 0; m <= m_max; ++m) {

		f128 df1 = double_factorial_128( (i+j-2*m) + k - l - 1);
		f128 df2 = double_factorial_128( (i+j-2*m) - k + l - 1);
		f128 df3 = double_factorial_128(-(i+j-2*m) + k + l - 1);

		f128 b1 = n_choose_k(i,m);
		f128 b2 = n_choose_k(j,m);

		f128 fm = factorial_128(m);
		//f128 f2 = factorial_128(i-m);
		//f128 f3 = factorial_128(j-m);

		//sum += (pow(2,m)/f1)*(1/f2)*(1/f3)*(df1*df2*df3);
		sum += (b1*ni)*(b2*nj)*(fm*nk)*(pow(2,m)*nl)*(df1*df2*df3);
	}

	const f64 sqrt_omega_over_2pi = sqrt(OMEGA/(2.0*M_PI));
	int raised = fetestexcept(FE_ALL_EXCEPT);
	if (raised) {
		sbmf_log_error("FPE ERROR!!");
		feclearexcept(FE_ALL_EXCEPT);
	}
	return sqrt_omega_over_2pi * sum;

	//return sqrt_omega_over_2pi*sqrt(((fi/fk)*(fj/fl)))/pow(2,(i+j+k+l)/2)*sum;
}
