#pragma once

#include "common.h"
#include <string.h>
#include <stdlib.h>

typedef struct {
	i32 order;
	i32 size;
	i32 bandcount;
	c64* bands;
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
static inline bandmat generate_fd_matrix(i32 n, i32 dim, f64 ds[]) {
	i32 k = pow(2,dim);

	bandmat b;
	b.size = pow(n, dim);
	b.order = b.size*b.size;
	b.bandcount = pow(n,dim-1)+1;
	b.bands = malloc(b.bandcount*b.size*sizeof(c64));
	memset(b.bands, 0, b.bandcount*b.size*sizeof(c64));

	// Setup main diagonal
	for (i32 i = (b.bandcount-1)*b.size; i < b.bandcount*b.size; ++i) {
		b.bands[i] = -k/(ds[0]*ds[0]);
	}

	// Setup off-diagonal elements
	for (i32 i = 1; i <= dim; ++i) {
		i32 bandindex = pow(n, i-1);
		for (i32 j = 0; j < b.size; ++j) {
			b.bands[j + (b.bandcount-1 - bandindex)*b.size] = 1/(ds[i-1]*ds[i-1]);
		}
	}

	return b;
}

