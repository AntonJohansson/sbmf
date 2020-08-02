#pragma once

#include <sbmf/common/common.h>
#include <sbmf/common/matrix.h>
#include <sbmf/debug/log.h>
#include <assert.h>

#define FACTORIAL_MAX_N 20
static u64 factorial_lookup[21] = {
	[0]  = 1,
	[1]  = 1,
	[2]  = 2,
	[3]  = 6,
	[4]  = 24,
	[5]  = 120,
	[6]  = 720,
	[7]  = 5040,
	[8]  = 40320,
	[9]  = 362880,
	[10] = 3628800,
	[11] = 39916800,
	[12] = 479001600,
	[13] = 6227020800,
	[14] = 87178291200,
	[15] = 1307674368000,
	[16] = 20922789888000,
	[17] = 355687428096000,
	[18] = 6402373705728000,
	[19] = 121645100408832000,
	[20] = 2432902008176640000
};


static inline u32 szudzik_pairing(u32 a, u32 b){
	return a >= b ? a*a + a + b : b*b + b + a;
}

static inline f64 factorial(u32 n) {
	if (n <= FACTORIAL_MAX_N) {
		return factorial_lookup[n];
	} else {
		return sqrt(M_2_PI*n) * pow(n/exp(1), n);
	}
}

static inline f64 ho_eigenfunction(i32 states[], f64 point[], i32 dims) {
	f64 prod = 1.0;
	static const f64 pi_factor = 1.0/pow(M_PI,0.25);
	for (i32 i = 0; i < dims; ++i) {
		i32 n = states[i];
		assert(n < 270);

		f64 x = point[i];

		f64 factorial_value = factorial(n);
		f64 normalization_factor = exp(-x*x/2.0) * pi_factor / sqrt(pow(2,n) * factorial_value);

		{
			f64 H0 = 1;
			f64 H1 = 2*x;
			f64 HN = 0.0;

			if (n == 0) {
				HN = H0;
			} else if (n == 1) {
				HN = H1;
			} else {
				for (i32 i = 2; i <= n; ++i) {
					HN = 2*x*H1 - 2*(i-1)*H0;
					H0 = H1;
					H1 = HN;
				}
			}

			prod *= normalization_factor*HN;
		}
	}
	return prod;
}

static f64 ho_eigenfunction_new(i32 states[], f64 point[], i32 dims) {
	f64 prod = 1.0;
	static const f64 pi_factor = 1.0/pow(M_PI,0.25);
	for (i32 i = 0; i < dims; ++i) {
		i32 n = states[i];
		assert(n < 270);

		f64 x = point[i];

		f64 H[3] = {1,1,1};
		for (i32 i = 1; i <= n; ++i) {
			H[2] = 2*(x*H[1] - (i-1)*H[0]);
			H[0] = H[1];
			H[1] = H[2];
		}

		f64 normalization_factor = exp(-x*x/2.0) * pi_factor / sqrt(pow(2,n) * factorial(n));
		prod *= normalization_factor*H[2];
	}
	return prod;
}

static void ho_eigenfunction_new_vec(u32 bases_per_dim[], f64 point[], u32 dims, f64* out) {
	static const f64 pi_factor = 1.0/pow(M_PI,0.25);
//
//	{
//		//u32 n[dims];
//		//memset(n, 0, dims*sizeof(u32));
//
//		//u32 tot = 1;
//		//for (u32 i = 0; i < dims; ++i)
//		//	tot *= s1[i]-s0[i];
//
//		// (3,3,3)
//		//
//		// 0  - 0,0,0 - 1
//		// 1  - 0,0,1 - 3
//		// 2  - 0,0,2 - 3
//		// 3  - 0,1,1 - 3
//		// 4  - 0,1,2 - 6
//		// 5  - 0,2,2 - 3
//		// 6  - 1,1,1 - 1
//		// 7  - 1,1,2 - 3
//		// 8  - 1,2,2 - 3
//		// 9  - 2,2,2 - 1
//		//
//		// 0  - 0,0,0,0 - 1
//		// 1  - 0,0,0,1 - 4
//		// 2  - 0,0,0,2 - 4
//		// 3  - 0,0,1,1 - 6
//		// 4  - 0,0,1,2 - 24
//		// 5  - 0,0,2,2 - 6
//		// 6  - 0,1,1,1 - 4
//		// 7  - 0,1,1,2 - 6
//		// 8  - 0,1,2,2 - 6
//		// 9  - 0,2,2,2 - 4
//		//
//		// 10 - 1,1,1,1 - 1
//		// 11 - 1,1,1,2 - 4
//		// 12 - 1,1,2,2 - 6
//		// 13 - 1,2,2,2 - 4
//		// 14 - 2,2,2,2 - 1
//		//
//		// 3*3*3*3/6 = 27*27/6 = 13.5
//		//
//		// 0,1,2
//		// 1,0,2
//		// 1,2,0
//		//
//		// 2,1,0
//		// 2,0,1
//		// 0,2,1
//
//		u32 n[3] = {0};
//		i32 index = 3-1;
//		for (u32 i = 0; i <= 3*3*3/3; ++i) {
//			// k = (index / (1))   		% l1
//			// j = (index / (l1))  		% l2
//			// i = (index / (l1*l2))  	% l3
//			// ...
//
//			log_info("%u\t- (%u,%u,%u)", i ,n[0], n[1], n[2]);
//
//			while (n[index]++ >= 2) {
//				n[index--] = 0;
//			}
//			while (index < 3) {
//				n[index+1] = n[index];
//				index++;
//			}
//			index = 3-1;
//
//			//u32 prodlen = 1;
//			//for (i32 j = dims-1; j >= 0; --j) {
//			//	n[j] = (i/prodlen) % (s1[j]-s0[j]);
//			//	prodlen *= (s1[j]-s0[j]);
//			//}
//
//			//for (u32 j = 0; j < dims; ++j) {
//			//	log_info("n: %d", n[j]);
//			//}
//			//log_info("\n");
//
//			//for (i32 n = 0; n < g.dimensions; ++n) {
//			//	g.points[g.dimensions*index + n] = g.mins[n] + indices[n]*g.deltas[n];
//			//}
//		}
//	}
//

	u32 tot = 1;
	for (u32 i = 0; i < dims; ++i)
		tot *= bases_per_dim[i];

	u32 n[dims];
	memset(n, 0, dims*sizeof(u32));

	for (u32 i = 0; i < tot; ++i) {
		for (u32 j = 0; j < dims; ++j)
			log_info("%u", n[j]);
		log_info("");

		{
			u32 index = dims-1;
			while (n[index]++ >= bases_per_dim[index]-1) {
				n[index--] = 0;
			}
		}
	}

	// Combination of states u00 u01 u10 u11:
	// 		(u0*u0) + (u0*u1) + (u1*u0) + (u1*u1)

	// So this is what we want to do for a combination of states
	// {n0,n1,n2,...}, but we still need to loop over the n's.
	for (u32 i = 0; i < dims; ++i) {
		u32 n = bases_per_dim[i];

		// maximum supported by 64-bit float, otherwise
		// errorso occur in hermite poly computation.
		assert(n < 270);

		for (u32 j = 0; j < n; ++j) {
			out[j] = 1.0;
		}

		f64 x = point[i];
		f64 psi_factor = exp(-x*x/2.0)*pi_factor;

		f64 H[3] = {1,1,1};
		out[0] = psi_factor;

		for (u32 j = 1; j < n; ++j) {
			H[2] = 2*(x*H[1] - (j-1)*H[0]);
			H[0] = H[1];
			H[1] = H[2];
			out[j] *= psi_factor*H[2]/sqrt(pow(2,j) * factorial(j));
		}
	}
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
