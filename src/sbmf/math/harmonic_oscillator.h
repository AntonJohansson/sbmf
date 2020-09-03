#pragma once

#include <sbmf/sbmf.h>
#include "matrix.h"
#include <sbmf/debug/log.h>

#include <assert.h>

static inline f128 factorial_128(u32 n) {
	f128 prod = 1.0;
	f128 current_value = (f128) n;
	while (current_value > 0.0) {
		prod *= current_value;
		current_value -= 1.0;
	}
	return prod;
}

/* Computed the value of the harmonic oscillator eigenfunction
 * at the passed in point in space and state-space
 */
static f64 ho_eigenfunction(i32 states[], f64 point[], i32 dims) {
	f64 prod = 1.0;
	static const f64 pi_factor = 1.0/pow(M_PI,0.25);
	for (i32 i = 0; i < dims; ++i) {
		i32 n = states[i];
		assert(n < 270);

		f64 x = point[i];

		f64 H[3] = {1,1,1};
		for (i32 j = 1; j <= n; ++j) {
			H[2] = 2*(x*H[1] - (j-1)*H[0]);
			H[0] = H[1];
			H[1] = H[2];
		}

		f64 normalization_factor = exp(-x*x/2.0) * pi_factor / sqrt(pow(2,n) * factorial_128(n));
		prod *= normalization_factor*H[2];
	}
	return prod;
}

/* Same as above but computes the eigenfunction for multiple state-values. Might remove. */
static void ho_eigenfunction_vec(u32 bases_per_dim[], f64 point[], u32 dims, f64* out) {
	static const f64 pi_factor = 1.0/pow(M_PI,0.25);

	/* Compute total dimensionality and find maximum states */
	u32 tot = 1;
	u32 max_state = 0;
	for (u32 i = 0; i < dims; ++i) {
		tot *= bases_per_dim[i];
		if (bases_per_dim[i] > max_state)
			max_state = bases_per_dim[i];
	}

	/* Storage for cached results */
	u32 store_index = 0;
	f64 store[max_state];

	u32 n[dims];
	memset(n, 0, dims*sizeof(u32));

	/* Combination of states u00 u01 u10 u11:
	 * 		(u0*u0) + (u0*u1) + (u1*u0) + (u1*u1)
	 */

	f64 sum = 0.0;
	for (u32 i = 0; i < tot; ++i) {
		f64 prod = 1.0;
		/* Go over each n0,n1,n2,... and compute H_(n0)(x),... */
		for (u32 j = 0; j < dims; ++j) {

			f64 x = point[j];

			f64 H[3] = {1,1,1};
			if (store_index > n[j]) {
				H[2] = store[n[j]];
			} else {
				/* Compute hermite polynomial and store it */
				for (u32 k = 1; k <= n[j]; ++k) {
					H[2] = 2*(x*H[1] - (k-1)*H[0]);
					H[0] = H[1];
					H[1] = H[2];
				}
				store[n[j]] = H[2];
				store_index++;
			}

			f64 psi_factor = exp(-x*x/2.0)*pi_factor;
			prod *= psi_factor*H[2]/sqrt(pow(2,n[j]) * factorial_128(n[j]));
		}
		sum += prod;

		/* Computes indices n0,n1,n2,... from global index i */
		{
			if (i == tot-1)
				break;

			u32 index = dims-1;
			while (n[index]++ >= bases_per_dim[index]-1) {
				n[index--] = 0;
			}
		}
	}
	log_info("sum: %lf", sum);

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

static inline c64 hob_sample(c64* v, u32 n, f64 x) {
	c64 output = 0;
	for (u32 i = 0; i < n; ++i) {
		output += v[i] * ho_eigenfunction((i32[]){i}, &x, 1);
	}
	return output;
}

/* <m|H|n> = <m|T|n> + <m|V|n>
 * T = p^2/2
 * * := adjoint operator
 * p = i(a* - a)/sqrt(2)
 * p^2 = -(a*a* - a*a - aa* + aa)/2
 * <m|a*a*|n> = <m|a*sqrt(n+1)|n+1>
 * 			  = sqrt[(n+2)(n+1)]<m|n+2>
 * 			  = sqrt[(n+2)(n+1)] if m == n+2,
 * 			    0 otherwise
 *
 * <m|a*a|n> = <m|a*sqrt(n)|n-1>
 * 			 = n<m|n>
 * 			 = n if m == n, 0 otherwise
 *
 * <m|aa*|n> = <m|asqrt(n+1)|n+1>
 * 			 = (n+1)<m|n>
 * 			 = (n+1) if m == n, 0 otherwise
 *
 * <m|aa|n> = <m|asqrt(n)|n-1>
 * 			= sqrt[n(n-1)]<m|n-2>
 * 			= sqrt(n(n-1)) if m == n-2, 0 otherwise
 *
 * that is the T matrix will have the following layout
 *
 * 		x   x
 * 		  x   x
 * 		x   x   x
 * 		  x   x
 * 		    x   x
 *
 * and <m|V|n> will be calculated numerically.
 */

static inline hermitian_bandmat construct_ho_kinetic_matrix(u32 size) {
	/* We only really need 2 bands for this matrix:
	 * the main diagonal + one off-diagonal band, since
	 * it is symmetric.
	 *
	 * However, we might as well allocate the full upper triangle
	 * since we'll add <m|V|n> terms to the diagonals later on
	 * anyway.
	 */
	hermitian_bandmat T = {
		.base = mat_new(size,size),
		.bandcount = size,
		.size = size,
	};

	/* How do we compute the position in band storage?
	 *
	 *		(x,x) (x,x) (x,x) (0,3)		0	1	2	3
	 *		(x,x) (x,x) (0,2) (1,3)		4	5	6	7
	 *		(x,x) (0,1) (1,2) (2,3)		8	9	10	11
	 *		(0,0) (1,1) (2,2) (3,3)		12	13	14	15
	 */

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
