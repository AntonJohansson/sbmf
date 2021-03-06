struct eigen_result_real find_eigenpairs_full_real(struct symmetric_bandmat bm) {
	f64 offdiag_elements[(bm.size-1)];
	f64 colmaj_eigvecs[bm.size*bm.size];
	f64 temp_eigvals[bm.size];

	/* Needed cuz dsbtrd overwrites diagonal of its input */
	f64 data[bm.size*bm.size];
	memcpy(data, bm.data, bm.size*bm.size*sizeof(f64));

	struct eigen_result_real res = {
		.eigenvalues            = sbmf_stack_push(bm.size * sizeof(f64)),
		.eigenvectors           = sbmf_stack_push(bm.size * bm.size * sizeof(f64)),
		.num_eigenpairs         = bm.size,
		.points_per_eigenvector = bm.size,
	};

	/* Start by reducing symmetric bandmatrix to tri-diagonal mat. */
	{
		i32 err = LAPACKE_dsbtrd(LAPACK_ROW_MAJOR, 'V', 'U',
				bm.size,
				bm.bandcount-1,
				data,
				bm.size,
				temp_eigvals, offdiag_elements, colmaj_eigvecs,
				bm.size);

		assert(err == 0); // @TODO, handle these errors better
	}

	/* Solve eigenvalue problem via QR factorisation of tridiagonal matrix. */
	{
		i32 err = LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'V',
				bm.size,
				temp_eigvals, offdiag_elements, colmaj_eigvecs,
				bm.size);

		assert(err == 0); // @TODO, handle these errors better
	}

	/* Convert eigenvectors to row-major, as we expected them to be. */
	for (u32 r = 0; r < bm.size; ++r) {
		res.eigenvalues[r] = temp_eigvals[r];
		for (u32 c = 0; c < bm.size; ++c) {
			res.eigenvectors[r + c*bm.size] = colmaj_eigvecs[c + r*bm.size];
		}
	}

	return res;
}

/* Arpack helper functions to provide more helpful error messages. */
static const char* arpack_dsaupd_error_code_to_string(i32 err);
static const char* arpack_dseupd_error_code_to_string(i32 err);
static const char* arpack_which_eigenpairs_to_string(enum which_eigenpairs which);

/*
 * Acutal code to solve the sparse eigenvalue problem.
 */
struct eigen_result_real find_eigenpairs_sparse_real(struct symmetric_bandmat bm, u32 num_eigenvalues, enum which_eigenpairs which) {
	u32 memory_marker = sbmf_stack_marker();

	i32 ido = 0;
	i32 n = bm.size;
	const char* which_str;
	switch (which) {
		case EV_LARGEST_MAG:  which_str = "LM"; break;
		case EV_SMALLEST_MAG: which_str = "SM"; break;
		case EV_LARGEST:	  which_str = "LA"; break;
		case EV_SMALLEST:	  which_str = "SA"; break;
		case EV_BOTH:		  which_str = "BE"; break;
		default: {
			sbmf_log_error("find_eigenparis_sprase_real(...): Invalid eigenvalue selected");
			return (struct eigen_result_real){0};
		 }
	};

	i32 nev = num_eigenvalues; // number of eigenvalues to find
	if (2 + nev > n) {
		sbmf_log_error("eig_sparse_bandmat(...) failed:");
		sbmf_log_error("\t the number of requested eigenvalues nev=%d and the matrix N=%d size must satisfy:", nev, n);
		sbmf_log_error("\t\t 2 + nev <= N");

		sbmf_stack_free_to_marker(memory_marker);
		return (struct eigen_result_real){0};
	}

	/* relative error tol */
	f64 tol = 0.0;
	/* residual vectors */
	f64 resid[n];

	/*
	 * Needs to satisfy
	 * 		2+nev <= ncv <= n
	 */
	i32 ncv = (2+nev+n)/2;

	f64* v = (f64*)sbmf_stack_push(sizeof(f64)*n*ncv);
	i32 iparam[11] = {
		[0] = 1, 		/* Exact shift is used */
		[2] = 7000, 	/* Max iterations */
		[3] = 1, 		/* will not work without this for some reason */
		[6] = 1, 		/* Mode */
	};
	i32 ipntr[14];

	f64* workd = (f64*)sbmf_stack_push(sizeof(f64)*3*n);
	i32 lworkl = ncv*ncv + 8*ncv;
	f64* workl = (f64*)sbmf_stack_push(sizeof(f64)*lworkl);

	i32 info = 0;

	struct symmetric_bandmat input_mat = bm;

	while (ido != 99 && info == 0) {

		dsaupd_c(&ido, "I", n, which_str, nev, tol, resid, ncv, v, n, iparam, ipntr,
				workd, workl, lworkl, &info);

		if (info != 0) {
			sbmf_log_error("arpack dsaupd(...) error: (%d) %s", info, arpack_dsaupd_error_code_to_string(info));
			break;
		}

		if (ido == -1 || ido == 1) {
			symmetric_bandmat_mulv(&workd[ipntr[1]-1], input_mat, &workd[ipntr[0]-1]);
		}
	}

	/* Convergence or error */
	if (info == 0) {
		i32 select[n]; // not used
		f64* d  = (f64*)sbmf_stack_push(sizeof(f64)*(nev));
		f64* z  = (f64*)sbmf_stack_push(sizeof(f64)*n*(nev));
		f64 sigma = 0; // not used

		dseupd_c(true, "A", select, d, z, n, sigma,
				"I", n, which_str, nev, tol, resid, ncv, v, n, iparam, ipntr, workd, workl, lworkl, &info);

		if (info != 0) {
			sbmf_log_error("arpack dseupd(...) error: (%d) %s", info, arpack_dseupd_error_code_to_string(info));
		} else {
#if 0
			i32 nconv = iparam[4];
			sbmf_log_info("arpack dsaupd(...) convergence info:");
			sbmf_log_info("\t Matrix size: %d", n);
			sbmf_log_info("\t Number of requested eigenvalues: %d", nev);
			sbmf_log_info("\t Number of generated Arnoldi vectors: %d", n);
			sbmf_log_info("\t Number of converged Ritz values: %d", nconv);
			sbmf_log_info("\t Sought eigenvalues: %s", which_str);
			sbmf_log_info("\t Number of iterations: %d", iparam[2]);
			sbmf_log_info("\t Number of OP*x: %d", iparam[8]);
			sbmf_log_info("\t Convergence tolerance: %lf", tol);
			sbmf_log_info("\t Eigenvalue residuals:");
			for (i32 i = 0; i < nconv; ++i) {
				c64 ax[n];
				complex_symmetric_bandmat_mulv(ax, input_mat, &z[i*n]);
				c64 neg_eigenvalue = -d[i];
				cblas_zaxpy(n, &neg_eigenvalue, &z[i*n], 1, ax, 1);
				sbmf_log_info("\t\t %lf + %lfi -- %e", creal(d[i]), cimag(d[i]), cblas_dznrm2(n, ax, 1));
			}
#endif

			//sbmf_stack_free_to_marker(memory_marker);

			struct eigen_result_real res = {
				.eigenvalues  = (f64*)sbmf_stack_push(sizeof(f64)*(nev)),
				.eigenvectors = (f64*)sbmf_stack_push(sizeof(f64)*bm.size*(nev)),
				.num_eigenpairs = nev,
				.points_per_eigenvector = bm.size,
			};
			memcpy(res.eigenvalues,  d, sizeof(f64)*nev);
			memcpy(res.eigenvectors, z, sizeof(f64)*bm.size*nev);

			return res;
		}
	}

	sbmf_log_error("not sure when this is called");

	sbmf_stack_free_to_marker(memory_marker);

	return (struct eigen_result_real){0};
}

/* See dsaupd.f for more info */
static inline const char* arpack_dsaupd_error_code_to_string(i32 err) {
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

/* See the dseupd.f for more info */
static inline const char* arpack_dseupd_error_code_to_string(i32 err) {
	switch (err) {
		case 0:  return "Normal exit.";
		case 1:  return "The Schur form computed by LAPACK routine csheqr could not be reordered by LAPACK routine ztrsen.  Re-enter subroutine dseupd with IPARAM(5)=NCV and increase the size of the array D to have dimension at least dimension NCV and allocate at least NCV columns for Z. NOTE: Not necessary if Z and V share the same space. Please notify the authors if this error occurs.";
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

static inline const char* arpack_which_eigenpairs_to_string(enum which_eigenpairs which) {
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
