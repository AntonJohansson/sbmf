#include "eigenproblem.h"

#include <lapacke.h>
#include <cblas.h>
#include <arpack/arpack.h>

#include <assert.h>
#include <stdlib.h> // malloc
#include <string.h> // memcpy, memset
#include <stdio.h> // printf


static inline void mat_transpose(c64* ans, c64* m, i32 rows, i32 cols) {
	c64 temp[rows*cols];
	memcpy(temp, m, sizeof(c64)*rows*cols);
	for (i32 c = 0; c < cols; ++c) {
		for (i32 r = 0; r < rows; ++r) {
			temp[r + c*rows] = m[c + r*cols];
		}
	}
	memcpy(ans, temp, sizeof(c64)*rows*cols);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void eig_dense_symetric_upper_tridiag_bandmat(bandmat bm, f64* out_eigvals, c64* out_eigvecs) {
	f64* offdiag_elements 	= malloc((bm.size - 1)*sizeof(f64));
	c64* colmaj_eigvecs = malloc((bm.size*bm.size)*sizeof(c64));

	i32 res = 0;

	// Start by reducing (complex hermitian) bandmatrix to tri-diagonal mat.
	{
		res = LAPACKE_zhbtrd(LAPACK_ROW_MAJOR, 'V', 'U', 
				bm.size,
				bm.super_diags, 
				bm.bands, 
				bm.size, 
				out_eigvals, offdiag_elements, colmaj_eigvecs, 
				bm.size);
		
		assert(res == 0); // @TODO, handle these errors better
	}

	// Solve eigenvalue problem via QR factorisation of tridiagonal matrix
	{
		res = LAPACKE_zsteqr(LAPACK_ROW_MAJOR, 'V', 
				bm.size, 
				out_eigvals, offdiag_elements, colmaj_eigvecs, 
				bm.size);

		assert(res == 0); // @TODO, handle these errors better
	}

	// Convert eigenvectors to row-major, as we expected them to be.
	mat_transpose(out_eigvecs, colmaj_eigvecs, bm.size, bm.size);
	
	free(offdiag_elements);
	free(colmaj_eigvecs);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static const char* arpack_znaupd_error_code_to_string(i32 err);
static const char* arpack_zneupd_error_code_to_string(i32 err);
static const char* which_eigenpairs_to_arpack_string(which_eigenpairs which);

void eig_sparse_bandmat(bandmat bm, u32 num_eigenvalues, which_eigenpairs which_pairs, c64* out_eigvals, c64* out_eigvecs) {
	// BEGIN PARAMS FOR ZNAUPD
	i32 ido = 0;
	char* bmat = "I";
	i32 n = bm.size;//*bm.size; // order of matrix
	const char* which = which_eigenpairs_to_arpack_string(which_pairs);
	i32 nev = num_eigenvalues; // number of eigenvalues to find
	f64 tol = 0.0; // relative error tolerance to achieve.
	c64 resid[n]; // residual vectors
	memset(resid, 0, sizeof(c64)*n);
	i32 ncv = n; // number of columns in V
	i32 ldv = n;
	c64 v[ncv*ldv];
	i32 iparam[11] = {
		[0] = 1, 		// Exact shift is used
		[2] = 300, 	// Max iterations
		[3] = 1, 		// will not work without this for some reason
		[6] = 1, 		// Mode
	};
	i32 ipntr[14];

	c64 workd[3*n];
	memset(workd, 0, sizeof(c64)*3*n);
	
	i32 lworkl = (3*ncv)*(3*ncv) + 5*ncv;
	c64 workl[lworkl];
	memset(workl, 0, sizeof(c64)*lworkl);

	f64 rwork[ncv];
	memset(rwork, 0, sizeof(f64)*ncv);
	i32 info = 0;
	// END PARAMS FOR ZNAUPD

	// BEGIN PARAMS FOR ZNEUPD
	bool rvec = true;
	char* howmny = "A";
	i32 select[ncv]; // not used
	c64 d[nev+1];
	c64 z[n*nev];
	i32 ldz = n;
	c64 sigma = 0; // not used
	c64 workev[2*ncv];
	// END PARAMS FOR ZNEUPD

	i32 kl = bm.super_diags; // number of superdiags
	i32 ku = bm.sub_diags; // number of subdiags
	i32 lda = kl + ku + 1;
			
	c64 one = 1, zero = 0;

	c64 a[lda*n];
	c64 m[lda*n];
	// Convert bands to col major ans save in a
	mat_transpose(a, bm.bands, lda, n);
 
	while (ido != 99 && info == 0) {
		znaupd_c(&ido, bmat, n, which, nev, tol, resid, ncv,
				v, ldv, iparam, ipntr, workd, workl, lworkl,
				rwork, &info);
		//printf("ido == %d\n", ido);

		if (ido == -1 || ido == 1) {
			cblas_zgbmv(CblasColMajor, CblasNoTrans, n, n, kl, ku, &one, a, 
					lda, &workd[ipntr[0]-1], 1, &zero, 
					&workd[ipntr[1]-1], 1);
		} else if (ido == 2) {
			cblas_zgbmv(CblasColMajor, CblasNoTrans, n, n, kl, ku, &one, m, 
					lda, &workd[ipntr[0]-1], 1, &zero, 
					&workd[ipntr[1]-1], 1);
		} else {
			// Convergence or errro
			if (info != 0) {
				// Error, just print the code
				printf("arpack znaupd error: %d (%s)\n", info, arpack_znaupd_error_code_to_string(info));
			} else {
				zneupd_c(rvec, howmny, select, d, z, ldz, sigma,
						workev, bmat, n, which, nev, tol,
						resid, ncv, v, ldv, iparam, ipntr, workd,
						workl, lworkl, rwork, &info);
				if (info != 0) {
					// Error just print code
					printf("arpack zneupd error: %d (%s)\n", info, arpack_zneupd_error_code_to_string(info));
				}
			}
		}
	}

	memcpy(out_eigvals, d, sizeof(c64)*(nev+1));
	mat_transpose(out_eigvecs, z, n, nev);
}

static inline const char* arpack_znaupd_error_code_to_string(i32 err) {
	switch (err) {
		case 0: return "Normal exit.";
		case 1: return "Maximum number of iterations taken.  All possible eigenvalues of OP has been found. IPARAM(5)  returns the number of wanted converged Ritz values.  ";
		case 2: return "No longer an informational error. Deprecated starting with release 2 of ARPACK.";
		case 3: return "No shifts could be applied during a cycle of the Implicitly restarted Arnoldi iteration. One possibility is to increase the size of NCV relative to NEV.  See remark 4 below.";
		case -1: return "N must be positive.";
		case -2: return "NEV must be positive.";
		case -3: return "NCV-NEV >case 2 and less than or equal to N.";
		case -4: return "The maximum number of Arnoldi update iteration must be greater than zero.";
		case -5: return "WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'";
		case -6: return "BMAT must be one of 'I' or 'G'.";
		case -7: return "Length of private work array is not sufficient.";
		case -8: return "Error return from LAPACK eigenvalue calculation;";
		case -9: return "Starting vector is zero.";
		case -10: return "IPARAM(7) must be 1,2,3.";
		case -11: return "IPARAM(7) case 1 and BMAT case 'G' are incompatable.";
		case -12: return "IPARAM(1) must be equal to 0 or 1.";
		case -9999: return "Could not build an Arnoldi factorization.  User input error highly likely.  Please check actual array dimensions and layout.  IPARAM(5) returns the size of the current Arnoldi factorization.";
		default: return "Unknown error code";	
	}
}
static inline const char* arpack_zneupd_error_code_to_string(i32 err) {
	switch (err) {
		case 0:  return "Normal exit.";
		case 1:  return "The Schur form computed by LAPACK routine csheqr could not be reordered by LAPACK routine ztrsen.  Re-enter subroutine zneupd with IPARAM(5)=NCV and increase the size of the array D to have dimension at least dimension NCV and allocate at least NCV columns for Z. NOTE: Not necessary if Z and V share the same space. Please notify the authors if this error occurs.";
		case -1: return "N must be positive.";
		case -2: return "NEV must be positive.";
		case -3: return "NCV-NEV >= 2 and less than or equal to N.";
		case -5: return "WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'";
		case -6: return "BMAT must be one of 'I' or 'G'.";
		case -7: return "Length of private work WORKL array is not sufficient.";
		case -8: return "Error return from LAPACK eigenvalue calculation. This should never happened.";
		case -9: return "Error return from calculation of eigenvectors. Informational error from LAPACK routine ztrevc.";
		case -10: return "IPARAM(7) must be 1,2,3";
		case -11: return "IPARAM(7) = 1 and BMAT = 'G' are incompatible.";
		case -12: return "HOWMNY = 'S' not yet implemented";
		case -13: return "HOWMNY must be one of 'A' or 'P' if RVEC = .true.";
		case -14: return "ZNAUPD did not find any eigenvalues to sufficient accuracy.";
		default: return "Unknown error code";	
	}
}
static inline const char* which_eigenpairs_to_arpack_string(which_eigenpairs which) {
	switch(which) {
		case EV_LARGEST_MAG:		return "LM";
		case EV_SMALLEST_MAG:		return "SM";
		case EV_LARGEST_RE:			return "LR";
		case EV_SMALLEST_RE:		return "SR";
		case EV_LARGEST_IM:			return "LI";
		case EV_SMALLEST_IM:		return "SI";
		default: 								return "unknown (which)";
	};
}
