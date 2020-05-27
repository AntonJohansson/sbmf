#include "eigenproblem.h"

#include <lapacke.h>
#include <cblas.h>
#include <arpack/arpack.h>

#include <assert.h>
#include <stdlib.h> // malloc
#include <string.h> // memcpy, memset
#include <stdio.h> // printf

#include <sbmf/memory/stack_allocator.h>

////////////////////////////////////////////////////////////////////////////////////////////////////

void eig_dense_symetric_upper_tridiag_bandmat(hermitian_bandmat bm, f64* out_eigvals, c64* out_eigvecs) {
	f64* offdiag_elements = (f64*)sa_push(_sbmf.main_stack, (bm.size-1)*sizeof(f64));
	c64* colmaj_eigvecs = (c64*)sa_push(_sbmf.main_stack, bm.size*bm.size*sizeof(c64));

	i32 res = 0;

	// Start by reducing (complex hermitian) bandmatrix to tri-diagonal mat.
	{
		res = LAPACKE_zhbtrd(LAPACK_ROW_MAJOR, 'V', 'U', 
				bm.size,
				bm.bandcount-1, 
				bm.base.data, 
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
	//TODO mat_transpose(out_eigvecs, colmaj_eigvecs, bm.size, bm.size);
	for (i32 r = 0; r < bm.size; ++r) {
		for (i32 c = 0; c < bm.size; ++c) {
			out_eigvecs[r + c*bm.size] = colmaj_eigvecs[c + r*bm.size];
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static const char* arpack_znaupd_error_code_to_string(i32 err);
static const char* arpack_zneupd_error_code_to_string(i32 err);
static const char* which_eigenpairs_to_arpack_string(which_eigenpairs which);

void eig_sparse_bandmat(hermitian_bandmat bm, u32 num_eigenvalues, which_eigenpairs which_pairs, c64* out_eigvals, c64* out_eigvecs) {
	i32 ido = 0;
	i32 n = bm.size;
	const char* which = which_eigenpairs_to_arpack_string(which_pairs);
	i32 nev = num_eigenvalues; // number of eigenvalues to find
	f64 tol = 0.1; // relative error tolerance to achieve.
	c64 resid[n]; // residual vectors

	i32 ncv = nev + 2;
	c64* v = (c64*)sa_push(_sbmf.main_stack, sizeof(c64)*n*ncv);
	i32 iparam[11] = {
		[0] = 1, 		// Exact shift is used
		[2] = 1000, 	// Max iterations
		[3] = 1, 		// will not work without this for some reason
		[6] = 1, 		// Mode
	};
	i32 ipntr[14];

	c64* workd = (c64*)sa_push(_sbmf.main_stack, sizeof(c64)*3*n);
	i32 lworkl = 3*ncv*ncv + 5*ncv;
	c64* workl = (c64*)sa_push(_sbmf.main_stack, sizeof(c64)*lworkl);
	f64 rwork[ncv];

	i32 info = 0;

	//c64 one = 1, zero = 0;
	//i32 num_of_diags = bm.super_diags + bm.sub_diags + 1; // 1 for main diag
	c64* data = (c64*)sa_push(_sbmf.main_stack, sizeof(c64)*bm.bandcount*n);
	hermitian_bandmat input_mat = {
		.base = {
			.rows = bm.bandcount,
			.cols = n,
			.data = data,
		},
		.bandcount = bm.bandcount,
		.size = bm.size,
	};
	mat_transpose(&input_mat.base, bm.base);
 
	while (ido != 99 && info == 0) {
		znaupd_c(&ido, "I", n, which, nev, tol, resid, ncv, v, n, iparam, ipntr, 
				workd, workl, lworkl, rwork, &info);

		if (info != 0) {
			printf("arpack znaupd error: %d (%s)\n", info, arpack_znaupd_error_code_to_string(info));
			break;
		}

		if (ido == -1 || ido == 1) {
			complex_hermitian_bandmat_mulv(&workd[ipntr[1]-1], input_mat, &workd[ipntr[0]-1]);
		}
	}

	// Convergence or error
	if (info == 0) {
		i32 select[n]; // not used
		//c64 d[nev+1];
		//c64 z[n*nev];
		c64* d = (c64*)sa_push(_sbmf.main_stack, sizeof(c64)*(nev+1));
		c64* z = (c64*)sa_push(_sbmf.main_stack, sizeof(c64)*n*nev);
		c64 sigma = 0; // not used
		//c64 workev[2*n];
		c64* workev = (c64*)sa_push(_sbmf.main_stack, sizeof(c64)*2*n);

		zneupd_c(true, "A", select, d, z, n, sigma, workev, 
				"I", n, which, nev, tol, resid, ncv, v, n, iparam, ipntr, workd, workl, lworkl, rwork, &info);

		if (info != 0) {
			printf("arpack zneupd error: %d (%s)\n", info, arpack_zneupd_error_code_to_string(info));
		} else {
			// Print info
			i32 nconv = iparam[4];
			printf("\nConvergence info:\n");
			printf("\t Matrix size: %d\n", n);
			printf("\t Number of requested eigenvalues: %d\n", nev);
			printf("\t Number of generated Arnoldi vectors: %d\n", n);
			printf("\t Number of converged Ritz values: %d\n", nconv);
			printf("\t Sought eigenvalues: %s\n", which);
			printf("\t Number of iterations: %d\n", iparam[2]);
			printf("\t Number of OP*x: %d\n", iparam[8]);
			printf("\t Convergence tolerance: %lf\n", tol);

			printf("\nEigenvalue -- Residual\n");
			for (i32 i = 0; i < nconv; ++i) {
				c64 ax[n];
				complex_hermitian_bandmat_mulv(ax, input_mat, &z[i*n]);
				c64 neg_eigenvalue = -d[i];
				cblas_zaxpy(n, &neg_eigenvalue, &z[i*n], 1, ax, 1);
				printf("\t %lf + %lfi -- %e\n", creal(d[i]), cimag(d[i]), cblas_dznrm2(n, ax, 1));
			}

			memcpy(out_eigvals, d, sizeof(c64)*(nev));
			mat_transpose_raw(out_eigvecs, z, n, nev);
		}
	}
}

static inline const char* arpack_znaupd_error_code_to_string(i32 err) {
	switch (err) {
		case 0: return "Normal exit.";
		case 1: return "Maximum number of iterations taken.  All possible eigenvalues of OP has been found. IPARAM(5)  returns the number of wanted converged Ritz values.  ";
		case 2: return "No longer an informational error. Deprecated starting with release 2 of ARPACK.";
		case 3: return "No shifts could be applied during a cycle of the Implicitly restarted Arnoldi iteration. One possibility is to increase the size of NCV relative to NEV.  See remark 4 below.";
		case -1: return "N must be positive.";
		case -2: return "NEV must be positive.";
		case -3: return "NCV-NEV >= 2 and less than or equal to N.";
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
