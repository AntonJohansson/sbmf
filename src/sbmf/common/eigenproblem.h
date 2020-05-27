#pragma once

#include "common.h"
#include "matrix.h"

typedef enum {
	EV_LARGEST_MAG 		= 0,
	EV_SMALLEST_MAG 	= 1,
	EV_LARGEST_RE			= 2,
	EV_SMALLEST_RE 		= 3,
	EV_LARGEST_IM 		= 4,
	EV_SMALLEST_IM		= 5
} which_eigenpairs;

extern void eig_dense_symetric_upper_tridiag_bandmat(hermitian_bandmat bm, f64* out_eigvals, c64* out_eigvecs);
extern void eig_sparse_bandmat(hermitian_bandmat bm, u32 num_eigenvalues, which_eigenpairs which_pairs, c64* out_eigvals, c64* out_eigvecs);
