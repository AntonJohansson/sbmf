#pragma once

#include <sbmf/common/common.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
#define FACTORIAL_MAX_N 20

static u64 factorial(u32 n) {
	assert(n <= FACTORIAL_MAX_N); // largest factorial supported by u64

	u64 res = 1;
	while (n > 1)
		res *= n--;

	return res;
}

static inline f64 hermite_poly(i32 n, f64 x) {
	f64 H0 = 1;
	f64 H1 = 2*x;

	if (n == 0)
		return H0;
	if (n == 1)
		return H1;

	f64 HN = 0.0;
	for (i32 i = 2; i <= n; ++i) {
		HN = 2*x*H1 - 2*(i-1)*H0;
		H0 = H1;
		H1 = HN;
	}

	return HN;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static inline f64 ho_eigenfunction(i32 states[], f64 point[], i32 dims) {
	f64 prod = 1.0;
	for (i32 i = 0; i < dims; ++i) {
		f64 normalization_factor = 0.0;
		if (states[i] <= FACTORIAL_MAX_N) {
			normalization_factor = 1.0/(pow(M_PI,0.25)*sqrt(pow(2,states[i])*factorial(states[i])));
		} else {
			f64 fac = sqrt(2*M_PI*states[i]) * pow( states[i]/exp(1), states[i]);
			normalization_factor = 1.0/(pow(M_PI,0.25)*sqrt(pow(2,states[i])*fac));
		}
		prod *= normalization_factor*exp(-point[i]*point[i]/2.0) * hermite_poly(states[i], point[i]);
	}
	return prod;
}

static inline f64 ho_eigenvalue(i32 states[], i32 dims) {
	f64 sum = 0.0;
	for (i32 i = 0; i < dims; ++i)
		sum += states[i];

	return sum + dims/2.0;
}

static inline f64 ho_potential(f64* v, i32 n, c64 u) {
	f64 temp = 0.0;
	for (i32 i = 0; i < n; ++i) {
		temp += v[i]*v[i];
	}
	return temp*0.5;
}
