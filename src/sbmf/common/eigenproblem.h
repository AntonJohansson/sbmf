#pragma once

#include "common.h"
#include "matrix.h"

typedef enum {
	EV_LARGEST_MAG 		= 0,
	EV_SMALLEST_MAG 	= 1,
	EV_LARGEST_RE		= 2,
	EV_SMALLEST_RE 		= 3,
	EV_LARGEST_IM 		= 4,
	EV_SMALLEST_IM		= 5
} which_eigenpairs;

typedef struct {
	c64* eigenvalues;
	c64* eigenvectors;
	u32 num_eigenpairs;
	u32 points_per_eigenvector;
} eig_result;

extern void eig_dense_symetric_upper_tridiag_bandmat(hermitian_bandmat bm, f64* out_eigvals, c64* out_eigvecs);
extern eig_result eig_sparse_bandmat(hermitian_bandmat bm, u32 num_eigenvalues, which_eigenpairs which_pairs);
