#pragma once

#include "matrix.h"
#include <sbmf/types.h>

enum which_eigenpairs {
	EV_LARGEST_MAG 		= 0,
	EV_SMALLEST_MAG 	= 1,
	EV_LARGEST_RE		= 2,
	EV_SMALLEST_RE 		= 3,
	EV_LARGEST_IM 		= 4,
	EV_SMALLEST_IM		= 5,
	EV_LARGEST			= 6,
	EV_SMALLEST 		= 7,
	EV_BOTH				= 8,
};

struct eigen_result {
	c64* eigenvalues;
	c64* eigenvectors;
	u32 num_eigenpairs;
	u32 points_per_eigenvector;
};

struct eigen_result_real {
	f64* eigenvalues;
	f64* eigenvectors;
	u32 num_eigenpairs;
	u32 points_per_eigenvector;
};

/* Find _all_ eigenpairs for a dense, symmetric,
 * upper tridiagonal matrix.
 */
struct eigen_result find_eigenpairs_full(struct complex_hermitian_bandmat bm);
struct eigen_result_real find_eigenpairs_full_real(struct hermitian_bandmat bm);

/* Find _some_ eigenpairs (specified by the enum which_eigenpairs)
 * for a dense, upper tridiagonal matrix.
 */
struct eigen_result_real find_eigenpairs_sparse_real(struct hermitian_bandmat bm, u32 num_eigenvalues, enum which_eigenpairs which);
struct eigen_result find_eigenpairs_sparse(struct complex_hermitian_bandmat bm, u32 num_eigenvalues, enum which_eigenpairs which);
