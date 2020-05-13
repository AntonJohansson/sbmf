#pragma once

#include <sbmf/common/common.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

static u64 factorial(u64 n) {
	assert(n <= 20); // largest factorial supported by u64

	if (n < 1) {
		return 1;
	} else {
		return n*factorial(n-1);
	}
}

static inline f64 hermite_poly(i32 n, f64 x) {
	// H_n(z) = (-1)^n * exp(z^2) * (d^n/dz^n) (exp(-z^2))
	//
	// Explicit form:
	// 		H_n(x) = n! sum_(m=0)^(floor(n/2)) (-1)^m/(m!*(n-2m)!) * (2x)^(n-2m)

	f64 sum = 0.0;
	i32 upper_bound = n/2;
	for (i32 m = 0; m <= upper_bound; ++m) {
		sum += pow(-1,m) / (factorial(m) * factorial(n-2*m)) * pow(2*x, n-2*m);
	}

	return factorial(n)*sum;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static inline f64 ho_eigenfunction(i32 states[], f64 point[], i32 dims) {
	f64 sum = 0.0;
	for (i32 i = 0; i < dims; ++i)
		sum += exp(-point[i]*point[i]/2.0) * hermite_poly(states[i], point[i]);
	return sum;
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
