#include "eigenproblem.h"

#include <lapacke.h>
#include <assert.h>
#include <stdlib.h>

void eigvp_bandmat(f64* out_eigval, c64* out_eigvec, bandmat bm) {
	f64* offdiag_elements 	= malloc((bm.size - 1)*sizeof(f64));
	c64* colmaj_eigvecs = malloc(bm.order*sizeof(c64));

	int res = 0;


	// Start by reducing (complex hermitian) bandmatrix to tri-diagonal mat.
	{
		res = LAPACKE_zhbtrd(LAPACK_ROW_MAJOR, 'V', 'U', 
				bm.size,
				bm.bandcount-1, 
				bm.bands, 
				bm.size, 
				out_eigval, offdiag_elements, colmaj_eigvecs, 
				bm.size);
		
		assert(res == 0);
	}

	// Solve eigenvalue problem via QR factorisation of tridiagonal matrix
	{
		res = LAPACKE_zsteqr(LAPACK_ROW_MAJOR, 'V', 
				bm.size, 
				out_eigval, offdiag_elements, colmaj_eigvecs, 
				bm.size);

		assert(res == 0);
	}

	// Convert eigenvectors to row-major, as we expected them to be.
	for (int i = 0; i < bm.size; ++i) {
		for (int j = 0; j < bm.size; ++j) {
			out_eigvec[i*bm.size + j] = colmaj_eigvecs[j*bm.size + i];
		}
	}
	
	free(offdiag_elements);
	free(colmaj_eigvecs);
}
