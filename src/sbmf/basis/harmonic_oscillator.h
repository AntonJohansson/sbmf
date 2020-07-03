#pragma once

#include <sbmf/common/common.h>
#include <sbmf/common/matrix.h>
#include <sbmf/debug/log.h>

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

// <m|H|n> = <m|T|n> + <m|V|n>
// T = p^2/2
// * := adjoint operator
// p = i(a* - a)/sqrt(2)
// p^2 = -(a*a* - a*a - aa* + aa)/2
// <m|a*a*|n> = <m|a*sqrt(n+1)|n+1>
// 			  = sqrt[(n+2)(n+1)]<m|n+2>
// 			  = sqrt[(n+2)(n+1)] if m == n+2,
// 			    0 otherwise
//
// <m|a*a|n> = <m|a*sqrt(n)|n-1>
// 			 = n<m|n>
// 			 = n if m == n, 0 otherwise
//
// <m|aa*|n> = <m|asqrt(n+1)|n+1>
// 			 = (n+1)<m|n>
// 			 = (n+1) if m == n, 0 otherwise
//
// <m|aa|n> = <m|asqrt(n)|n-1>
// 			= sqrt[n(n-1)]<m|n-2>
// 			= sqrt(n(n-1)) if m == n-2, 0 otherwise
//
// that is the T matrix will have the following layout
//
// 		x   x
// 		  x   x
// 		x   x   x
// 		  x   x
// 		    x   x
//
// and <m|V|n> will be calculated numerically.

static inline hermitian_bandmat construct_ho_kinetic_matrix(u32 size) {
	// We only really need 2 bands for this matrix:
	// the main diagonal + one off-diagonal band, since
	// it is symmetric.
	//
	// However, we might as well allocate the full upper triangle
	// since we'll add <m|V|n> terms to the diagonals later on
	// anyway.
	//
	hermitian_bandmat T = {
		.base = mat_new(size,size),
		.bandcount = size,
		.size = size,
	};

	// How do we compute the position in band storage?
	//
	//		(x,x) (x,x) (x,x) (0,3)		0	1	2	3
	//		(x,x) (x,x) (0,2) (1,3)		4	5	6	7
	//		(x,x) (0,1) (1,2) (2,3)		8	9	10	11
	//		(0,0) (1,1) (2,2) (3,3)		12	13	14	15

	for (i32 r = 0; r < size; ++r) {
		for (i32 c = r; c < size; ++c) {
			u32 i = (size-1)*(size-(c-r)) + r;
			if (r == c) {
				T.base.data[i] = (2*c + 1)/4.0;
			} else if (r == c - 2) {
				T.base.data[i] = -sqrt(c*(c-1))/4.0;
			} else {
				T.base.data[i] = 0;
			}
		}
	}

	return T;
}
