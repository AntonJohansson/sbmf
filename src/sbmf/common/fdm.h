#pragma once

#include "common.h"
#include <string.h>
#include <stdlib.h>

typedef struct {
	int_t order;
	int_t size;
	int_t bandcount;
	complex_t* bands;
} bandmat;

// Generates [n x n] finite difference matrix of order dim, with diagonal element -k.
// For example 
//
// 	generate_fd_matrix(5, 2, 1, (real_t[]){dx})
//
// gives the matrix
//
//		 	    [-2  1  0  0  0 ]
//		 	1   [ 1 -2  1  0  0 ]
//		 ---  [ 0  1 -2  1  0 ]
//		 dx^2 [ 0  0  1 -2  1 ]
//		 	    [ 0  0  0  1 -2 ]
//
// expressed as band array
//
// 	 bands = (1/dx^2) * [ 1  1  1  1  1 
//                       -2 -2 -2 -2 -2 ].
//
// For the 2D case
//
// 		generate_fd_matrix(3, 4, 2, (real_t[]){dx,dy}),
//
// gives the matrix
//
// 		[ -4  1  0  1  0             ]
// 		[  1 -4  1  0  1  0          ]
// 		[  0  1 -4  1  0  1  0       ]
// 		[  1  0  1 -4  1  0  1  0    ]
// 		[  0  1  0  1 -4  1  0  1  0 ]
// 		[     0  1  0  1 -4  1  0  1 ]
// 		[        0  1  0  1 -4  1  0 ]
// 		[           0  1  0  1 -4  1 ]
// 		[              0  1  0  1 -4 ]
//
// where the central tridiagonal is scaled by 1/dx^2 and the off-diagonal is scaled by
// 1/dy^2. The matrix is represented as the band array
//
//  	bands = [  1  1  1  1  1  1  1  1  1 
//  	           0  0  0  0  0  0  0  0  0
//							 1  1  1  1  1  1  1  1  1
//  	          -4 -4 -4 -4 -4 -4 -4 -4 -4 ],
//
// where it has been assumed that dx = dy = 1.
static inline bandmat generate_fd_matrix(int_t n, int_t k, int_t dim, real_t ds[]) {
	bandmat b;
	b.size = pow(n, dim);
	b.order = b.size*b.size;
	b.bandcount = pow(n,dim-1)+1;
	b.bands 	= malloc(b.bandcount*b.size*sizeof(complex_t));
	memset(b.bands, 0, b.bandcount*b.size*sizeof(complex_t));

	for (int_t i = (b.bandcount-1)*b.size; i < b.bandcount*b.size; ++i) {
		b.bands[i] = -k/(ds[0]*ds[0]);
	}

	for (int_t i = 1; i <= dim; ++i) {
		int_t bandindex = pow(n, i-1);
		for (int_t j = 0; j < b.size; ++j) {
			b.bands[j + (b.bandcount-1 - bandindex)*b.size] = 1/(ds[i-1]*ds[i-1]);
		}
	}

	return b;
}

