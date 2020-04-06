#pragma once

#include "common.h"
#include <stdio.h>

// Convenience macro to loop over a matrix (rows x cols) row-wise, i.e.
// in memory-order.
#define FOREACH_ROW(rows, cols, i,j) \
	for (int_t i = 0; i < rows; ++i) \
		for (int_t j = 0; j < cols; ++j)

// Convenience function to calculate index into the array corresponding to the
// (i,j):th element in a the matrix of size NxN.
static inline int_t mat_idx(int_t cols, int_t i, int_t j) { 
	return j + cols*i;
}

// Returns the element of maximum absolute value for the input matrix.
static inline real_t max_abs_complex(int_t rows, int_t cols, complex_t* m) {
	real_t max_norm = -DBL_MAX;
	FOREACH_ROW(rows,cols, i,j) {
		real_t d = cabs(m[mat_idx(rows,i,j)]);
		if (d > max_norm) {
			max_norm = d;
		}
	}
	return max_norm;
}

// Prints matrix of size (rows x cols) to stdout.
static inline void fprintcmat(FILE* fd, int_t rows, int_t cols, complex_t* m) {
	fprintf(fd, "matrix [%dx%d]\n", rows, cols);
	for (int_t i = 0; i < rows; ++i) {
		for (int_t j = 0; j < cols; ++j) {
			real_t cr = creal(m[mat_idx(cols,i,j)]);
			real_t ci = cimag(m[mat_idx(cols,i,j)]);
			fprintf(fd, "%lf%s%lfi", cr, (ci < 0.0) ? "-" : "+", fabs(ci));

			if (j < cols-1)
				fprintf(fd, "\t");
		}
		fprintf(fd, "\n");
	}
}

static inline void fprintrmat(FILE* fd, int_t rows, int_t cols, real_t* m) {
	fprintf(fd, "matrix [%dx%d]\n", rows, cols);
	for (int_t i = 0; i < rows; ++i) {
		for (int_t j = 0; j < cols; ++j) {
			fprintf(fd, "%lf", m[mat_idx(cols,i,j)]);
			if (j < cols-1)
				fprintf(fd, "\t");
		}
		fprintf(fd, "\n");
	}
}
