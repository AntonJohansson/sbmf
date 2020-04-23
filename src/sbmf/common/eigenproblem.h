#pragma once

#include "fdm.h" // for bandmat
#include <lapacke.h>
#include <assert.h>
#include <stdlib.h>

static inline void eigvp_bandmat(real_t* out_eigval, complex_t* out_eigvec, bandmat bm) {
	real_t* offdiag_elements = malloc((bm.size - 1)*sizeof(real_t));

	int res = 0;

	// Start by reducing (complex hermitian) bandmatrix to tri-diagonal mat.
	{
		res = LAPACKE_zhbtrd(LAPACK_ROW_MAJOR, 'V', 'U', 
				bm.size,
				bm.bandcount-1, 
				bm.bands, 
				bm.size, 
				out_eigval, offdiag_elements, out_eigvec, 
				bm.size);
		
		assert(res == 0);
	}

	// Solve eigenvalue problem via QR factorisation of tridiagonal matrix
	{
		res = LAPACKE_zsteqr(LAPACK_ROW_MAJOR, 'V', 
				bm.size, 
				out_eigval, offdiag_elements, out_eigvec, 
				bm.size);

		assert(res == 0);
	}

	free(offdiag_elements);
}

/*
	 {
	 real_t* outD = malloc(fdm.size*sizeof(real_t));
	 real_t* outE = malloc((fdm.size-1) * sizeof(real_t));
	 complex_t* outQ = malloc(fdm.order*sizeof(complex_t));

	 {
	 int res = 0;

// Reduce to tridiag-form
res = LAPACKE_zhbtrd(LAPACK_ROW_MAJOR, 'V', 'U', 
fdm.size,
fdm.bandcount-1, 
fdm.bands, 
fdm.size, 
outD, outE, outQ, 
fdm.size);

assert(res == 0);

// Solve eigenvalue problem
res = LAPACKE_zsteqr(LAPACK_ROW_MAJOR, 'V', 
fdm.size, 
outD, outE, outQ, 
fdm.size);

assert(res == 0);
}

// outD contains diagonal elements
// outE has been destroyed on exit
// outQ contains eigenvectors

float x[g.total_pointcount];
float u[g.total_pointcount];
for (int i = 1; i < 6; ++i) {
//normalize_wf(&outQ[i*g.total_pointcount], g.total_pointcount);
for (int j = 0; j < g.total_pointcount; ++j) {
x[j] = (float)g.points[j];
u[j] = (float)outD[i] +  5*(float)creal(outQ[i*g.total_pointcount + j]);
}
plt_1d(state, x,u, g.total_pointcount);
}
}
*/
