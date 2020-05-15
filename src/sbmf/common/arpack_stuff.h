#pragma once

static inline void eigvp_bandmat_(_bandmat bm, i32 number_of_eigenvalues, f64* out_eigvals, c64* out_eigvecs) {
////////////////////////////////////////////////////////////////////////////////////////////////////
		if (info == 0) {
			//nconv = iparam[5];
			//printf(
			//		"The size of the matrix is %d\n"
			//		"Number of eigenvalue requested is %d\n"
			//		"The number of Arnoldi vectors generated is %d\n"
			//		"The number of converged Ritz values is %d\n"
			//		"What portion of the spectrum %c%c\n"
			//		"The number of implicit Arnoldi update takes is %d\n"
			//		"The number of OP*x is %d\n"
			//		"The convergence tolerance is %lf\n",
			//		n, nev, ncv, which[0], which[1], iparam[3], iparam[9], tol);
		// Compute residual norm || Ax - lx ||
		//for (i32 j = 0; j < nconv; ++j) {
    //        cblas_zgbmv(CblasRowMajor, CblasNoTrans,  n, n, kl, ku, 
		//						&one,
		//						&a[kl+1,1], lda, 
		//						&v[1,j], 1, 
		//						&zero,
		//						&ax, 1);

		//				c64 alpha = -d[j];
		//				cblas_zaxpy(n, &alpha, &v[1,j], 1, &ax, 1);
    //        rd(j,1) = dble(d(j))
    //        rd(j,2) = dimag(d(j))
    //        rd(j,3) = dznrm2(n, ax, 1)
    //        rd(j,3) = rd(j,3) / dlapy2(rd(j,1),rd(j,2))
		//}
         //dmout(6, nconv, 3, rd, maxncv, -6,
         //         'Ritz values (Real,Imag) and relative residuals')
		} else {
			// Print znband error
		}
}
