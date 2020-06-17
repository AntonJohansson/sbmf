#define SNOW_ENABLED
#include "snow.h"

#include <sbmf/groundstate_solver/groundstate_solver.h>
#include <sbmf/common/profile.h>
#include <sbmf/common/fdm.h>
#include <sbmf/common/grid.h>

#include <sbmf/common/eigenproblem.h>

#include <plot/plot.h>

#include <sbmf/basis/harmonic_oscillator.h>

#define PLOT_SPARSE_1D_HO_FDM 0
#define PLOT_SPARSE_2D_HO_FDM 0

#define PLOT_SPARSE_1D_PB_FDM 0
#define PLOT_SPARSE_2D_PB_FDM 0

#define PLOT_SPARSE_1D_PB_HOBASIS 1

////////////////////////////////////////////////////////////////////////////////////////////////////
// TEST NUMERICAL INTEGRATION VIA QUADGK
#include <sbmf/quadgk.h>

static f64 integrand_cos(f64 x, void* data) {return cos(x);}
static f64 integrand_sin(f64 x, void* data) {return sin(x);}
static f64 integrand_x2(f64 x, void* data) { return x*x; }
static f64 integrand_expx(f64 x, void* data) { return exp(x); }
static f64 integrand_expnx(f64 x, void* data) { return exp(-x); }
static f64 integrand_expnabsx(f64 x, void* data) { return exp(-fabs(x)); }

describe(quadgk_integration) {
	integration_settings settings = {
		.order = 7,
		.abs_error_tol = 0,
		.rel_error_tol = 1e-4,
		.max_evals = 1e7,
	};

	integration_result res;

	it ("cosine 0 -> 1") {
		res = quadgk(integrand_cos, 0,1, settings);
		assert(res.integral - 0.84 < 0.01);
		assert(res.error < 1e-4);
	}

	it ("exp 0 -> 1") {
		res = quadgk(integrand_expx, 0,1, settings);
		assert(res.integral - 1.78 < 0.01);
		assert(res.error < 1e-4);
		assert(res.converged);
	}

	it ("sine 0 -> 2pi") {
		res = quadgk(integrand_sin, 0,2*M_PI, settings);
		assert(res.integral - 0.0 < 0.01);
		assert(res.error < 1e-4);
		assert(res.converged);
	}

	it ("x^2 -1 -> 1") {
		res = quadgk(integrand_x2, -1, 1, settings);
		assert(res.integral - 2.0/3.0 < 0.01);
		assert(res.error < 1e-4);
		assert(res.converged);
	}

	it ("exp(-x) 0 -> inf") {
		res = quadgk(integrand_expnx, 0, INFINITY, settings);
		//printf("\n%lf\n", res.integral);
		assert(res.integral - 1.0 < 0.01);
		assert(res.error < 1e-4);
		assert(res.converged);
	}

	it ("exp(-x) inf -> 0") {
		res = quadgk(integrand_expnx, INFINITY, 0, settings);
		assert(res.integral - (-1.0) < 0.01);
		assert(res.error < 1e-4);
		assert(res.converged);
	}

	it ("exp(x) -inf -> 0") {
		res = quadgk(integrand_expx, -INFINITY, 0, settings);
		assert(res.integral - (1.0) < 0.01);
		assert(res.error < 1e-4);
		assert(res.converged);
	}

	it ("exp(x) 0 -> -inf") {
		res = quadgk(integrand_expx, 0, -INFINITY, settings);
		//printf("\n%lf\n", res.integral);
		assert(res.integral - (-1.0) < 0.01);
		assert(res.error < 1e-4);
		assert(res.converged);
	}

	it ("exp(-|x|) -inf -> inf") {
		res = quadgk(integrand_expnabsx, -INFINITY, INFINITY, settings);
		assert(res.integral - (2.0) < 0.01);
		assert(res.error < 1e-4);
		assert(res.converged);
	}

	it ("exp(-|x|) inf -> -inf") {
		res = quadgk(integrand_expnabsx, -INFINITY, INFINITY, settings);
		//printf("\n%lf\n", res.integral);
		assert(res.integral - (-2.0) < 0.01);
		assert(res.error < 1e-4);
		assert(res.converged);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// TEST PRIORITY QUEUE
#include <sbmf/prioqueue.h>

static bool cmpint(void* a, void* b) {
	return (*(i32*)a > *(i32*)b);
}

describe(priority_queue) {
	it ("push items") {
		prioqueue* pq = prioqueue_new(10, sizeof(i32), cmpint);

		i32 i;
		i	= 0; prioqueue_push(pq, &i);
		i	= 1; prioqueue_push(pq, &i);
		i	= 2; prioqueue_push(pq, &i);
		i	= 3; prioqueue_push(pq, &i);
		i	= 4; prioqueue_push(pq, &i);
		i	= 5; prioqueue_push(pq, &i);
		i	= 6; prioqueue_push(pq, &i);
		i	= 7; prioqueue_push(pq, &i);
		i	= 8; prioqueue_push(pq, &i);
		i	= 9; prioqueue_push(pq, &i);

		{
			i32 expected[] = {9,8,7,6,5,4,3,2,1,0};
			for (i32 i = 0; i < ARRLEN(expected); ++i) {
				i32* ptr = (i32*)prioqueue_top(pq);
				//printf("%d\n", *ptr);
				asserteq(*ptr, expected[i]);
				prioqueue_pop(pq);
			}

			asserteq(pq->size, 0);
		}

		prioqueue_free(pq);
	}
	it ("push and pop items") {
		prioqueue* pq = prioqueue_new(10, sizeof(i32), cmpint);

		i32 i;
		i	= 0; prioqueue_push(pq, &i);
		i	= 1; prioqueue_push(pq, &i);
		i	= 2; prioqueue_push(pq, &i);
		i	= 3; prioqueue_push(pq, &i);
		i	= 4; prioqueue_push(pq, &i);
		i	= 5; prioqueue_push(pq, &i);
		i	= 6; prioqueue_push(pq, &i);
		i	= 7; prioqueue_push(pq, &i);
		i	= 8; prioqueue_push(pq, &i);
		i	= 9; prioqueue_push(pq, &i);

		prioqueue_pop(pq);
		prioqueue_pop(pq);
		prioqueue_pop(pq);

		i	= 10; prioqueue_push(pq, &i);
		i	= 12; prioqueue_push(pq, &i);
		i	= 11; prioqueue_push(pq, &i);

		{
			i32 expected[] = {12,11,10,6,5,4,3,2,1,0};
			for (i32 i = 0; i < ARRLEN(expected); ++i) {
				i32* ptr = (i32*)prioqueue_top(pq);
				//printf("%d\n", *ptr);
				asserteq(*ptr, expected[i]);
				prioqueue_pop(pq);
			}

			asserteq(pq->size, 0);
		}

		prioqueue_free(pq);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Hamiltonian solving via FDM

describe(fdm_ho) {
	it ("fdm construction 1D [5x5]") {
		sbmf_init();
		f64 deltas[5] = {1,1,1,1,1};
		hermitian_bandmat fdm = construct_finite_diff_mat(5, 1, deltas);
		c64 expected_ans[] = {
			 1, 1, 1, 1, 1,
			-2,-2,-2,-2,-2
		};
		for (i32 i = 0; i < fdm.base.rows*fdm.base.cols; ++i) {
			asserteq(fdm.base.data[i], expected_ans[i]);
		}
		sbmf_shutdown();
	}
	it ("1D -- full diag") {
		sbmf_init();

		const i32 N = 128;
		grid g = generate_grid(1,
				(f64[]){-5},
				(f64[]){ 5},
				(i32[]){ N});

		hermitian_bandmat fdm = construct_finite_diff_mat(N, g.dimensions, g.deltas);
		for (i32 i = 0; i < fdm.size*fdm.bandcount; ++i) {
			fdm.base.data[i] = -0.5*fdm.base.data[i];
		}
		for (i32 i = 0; i < fdm.size; ++i) {
			i32 idx = fdm.size*(fdm.bandcount-1) + i;
			fdm.base.data[idx] += (c64) ho_potential(&g.points[g.dimensions*i], g.dimensions, 0.0);
		}

		f64 out_eigvals[N];
		c64 out_eigvecs[N*N];
		eig_dense_symetric_upper_tridiag_bandmat(fdm, out_eigvals, out_eigvecs);

		sbmf_shutdown();
	}
	it ("1D -- sparse diag") {
		sbmf_init();
		const i32 N = 40;
		grid g = generate_grid(1,
				(f64[]){-5},
				(f64[]){ 5},
				(i32[]){ N});

		hermitian_bandmat fdm = construct_finite_diff_mat(N, g.dimensions, g.deltas);
		for (i32 i = 0; i < fdm.size*fdm.bandcount; ++i) {
			fdm.base.data[i] = -0.5*fdm.base.data[i];
		}
		for (i32 i = 0; i < fdm.size; ++i) {
			i32 idx = fdm.size*(fdm.bandcount-1) + i;
			fdm.base.data[idx] += (c64) ho_potential(&g.points[g.dimensions*i], g.dimensions, 0.0);
		}

		u32 eigenvalues_to_find = 3;
		eig_result res = eig_sparse_bandmat(fdm, eigenvalues_to_find, EV_SMALLEST_RE);

		// Plotting
		{
#if PLOT_SPARSE_1D_HO_FDM
			plot_init(800, 600, "1D HO FDM");
			sample_space lin = make_linspace(1,-5,5,N);
			f32 z[lin.pointcount];

			FOREACH_POINT_1D(lin, i)
				z[i] = ho_potential(&g.points[i], 1, 0);
			push_line_plot(&(plot_push_desc){
					.space = &lin,
					.data = z,
					.label = "potential",
					});

			for (u32 n = 0; n < res.num_eigenpairs; ++n) {
				FOREACH_POINT_1D(lin, i)
					z[i] = ho_eigenfunction((i32[]){n}, &g.points[i], 1);
				push_line_plot(&(plot_push_desc){
						.space = &lin,
						.data = z,
						.label = plot_snprintf("exact %d state", n),
						.offset = ho_eigenvalue((i32[]){n},1),
						});
			}

			for (u32 n = 0; n < res.num_eigenpairs; ++n) {
				FOREACH_POINT_1D(lin, i)
					z[i] = cabs(res.eigenvectors[n*fdm.size + i]);
				push_line_plot(&(plot_push_desc){
						.space = &lin,
						.data = z,
						.label = plot_snprintf("numerical %d state", n),
						.offset = res.eigenvalues[n],
						});
			}

			plot_update_until_closed();
			plot_shutdown();
#endif
		}

		sbmf_shutdown();
	}
	it ("2D") {
		sbmf_init();

		const i32 N = 32;
		grid g = generate_grid(2,
				(f64[]){-2.5,-2.5},
				(f64[]){ 2.5, 2.5},
				(i32[]){ N, N});

		hermitian_bandmat fdm = construct_finite_diff_mat(N, g.dimensions, g.deltas);
		for (i32 i = 0; i < fdm.size*fdm.bandcount; ++i) {
			fdm.base.data[i] = -0.5*fdm.base.data[i];
		}
		for (i32 i = 0; i < fdm.size; ++i) {
			i32 idx = fdm.size*(fdm.bandcount-1) + i;
			fdm.base.data[idx] += (c64) ho_potential(&g.points[g.dimensions*i], g.dimensions, 0.0);
		}

		u32 eigenvalues_to_find = 2;
		eig_result res = eig_sparse_bandmat(fdm, eigenvalues_to_find, EV_SMALLEST_RE);


		// Plotting
		{
#if PLOT_SPARSE_2D_HO_FDM
			plot_init(800, 600, "2D HO FDM");
			sample_space surf = make_linspace(g.dimensions, -2.5, 2.5, N);
			f32 z[surf.pointcount];

			FOREACH_POINT_1D(surf, i)
				z[i] = ho_potential(&g.points[g.dimensions*i], g.dimensions, 0.0);
			push_surface_plot(&(plot_push_desc){
					.space = &surf,
					.data = z,
					.label = "potential",
					});

			for (u32 n1 = 0; n1 < res.num_eigenpairs; ++n1) {
				for (u32 n2 = 0; n2 < res.num_eigenpairs; ++n2) {
					i32 states[] = {n1,n2};
					FOREACH_POINT_1D(surf, i)
						z[i] =  ho_eigenfunction(states, &g.points[g.dimensions*i], g.dimensions);
					push_surface_plot(&(plot_push_desc){
							.space = &surf,
							.data = z,
							.label = plot_snprintf("exact %d,%d state", n1,n2),
							.offset = ho_eigenvalue(states, g.dimensions),
							});
				}
			}

			for (u32 n = 0; n < res.num_eigenpairs; ++n) {
				FOREACH_POINT_1D(surf, i)
					z[i] = 100*creal(res.eigenvectors[n*fdm.size + i]);
				push_surface_plot(&(plot_push_desc){
						.space = &surf,
						.data = z,
						.label = plot_snprintf("numerical %d state", n),
						.offset = res.eigenvalues[n],
						});
			}

			plot_update_until_closed();
			plot_shutdown();
#endif
		}

		sbmf_shutdown();
	}
}

static f64 pb_potential(f64* point, u32 dims) {
	for (u32 i = 0; i < dims; ++i)
		if (abs(point[i]) >= 3.0)
			return 10.0;
	return 0.0;
}

describe(fdm_pb) {
	it ("1D") {
		sbmf_init();
		const i32 N = 40;
		grid g = generate_grid(1,
				(f64[]){-5},
				(f64[]){ 5},
				(i32[]){ N});

		hermitian_bandmat fdm = construct_finite_diff_mat(N, g.dimensions, g.deltas);
		for (i32 i = 0; i < fdm.size*fdm.bandcount; ++i) {
			fdm.base.data[i] = -0.5*fdm.base.data[i];
		}
		for (i32 i = 0; i < fdm.size; ++i) {
			i32 idx = fdm.size*(fdm.bandcount-1) + i;
			fdm.base.data[idx] += (c64) pb_potential(&g.points[g.dimensions*i], g.dimensions);
		}

		u32 eigenvalues_to_find = 3;
		eig_result res = eig_sparse_bandmat(fdm, eigenvalues_to_find, EV_SMALLEST_RE);

		// Plotting
		{
#if PLOT_SPARSE_1D_PB_FDM
			plot_init(800, 600, "1D PB FDM");

			sample_space sp = make_linspace(g.dimensions, -5, 5, N);
			f32 z[sp.pointcount];

			FOREACH_POINT_1D(sp, i)
				z[i] = pb_potential(&g.points[g.dimensions*i], g.dimensions);
			push_line_plot(&(plot_push_desc){
						.space = &sp,
						.data = z,
						.label = "potential",
					});

			for (u32 n = 0; n < res.num_eigenpairs; ++n) {
				FOREACH_POINT_1D(sp, i)
					z[i] = creal(res.eigenvectors[n*res.points_per_eigenvector + i]);
				push_line_plot(&(plot_push_desc){
						.space = &sp,
						.data = z,
						.label = plot_snprintf("numerical %d state", n),
						.offset = res.eigenvalues[n],
						});
			}

			plot_update_until_closed();
			plot_shutdown();
#endif
		}

		sbmf_shutdown();
	}
	it ("2D") {
		sbmf_init();
		const i32 N = 32;
		grid g = generate_grid(2,
				(f64[]){-5, -5},
				(f64[]){ 5,  5},
				(i32[]){ N,  N});

		hermitian_bandmat fdm = construct_finite_diff_mat(N, g.dimensions, g.deltas);
		for (i32 i = 0; i < fdm.size*fdm.bandcount; ++i) {
			fdm.base.data[i] = -0.5*fdm.base.data[i];
		}
		for (i32 i = 0; i < fdm.size; ++i) {
			i32 idx = fdm.size*(fdm.bandcount-1) + i;
			fdm.base.data[idx] += (c64) pb_potential(&g.points[g.dimensions*i], g.dimensions);
		}

		u32 eigenvalues_to_find = 4;
		eig_result res = eig_sparse_bandmat(fdm, eigenvalues_to_find, EV_SMALLEST_RE);

		// Plotting
		{
#if PLOT_SPARSE_2D_PB_FDM
			plot_init(800, 600, "2D PB FDM");

			sample_space sp = make_linspace(g.dimensions, -5, 5, N);
			f32 z[sp.pointcount];

			FOREACH_POINT_1D(sp, i)
				z[i] = pb_potential(&g.points[g.dimensions*i], g.dimensions);
			push_surface_plot(&(plot_push_desc){
					.space = &sp,
					.data = z,
					.label = "potential",
					});

			for (u32 n = 0; n < res.num_eigenpairs; ++n) {
				FOREACH_POINT_1D(sp, i)
					z[i] = 100*creal(res.eigenvectors[n*res.points_per_eigenvector + i]);
				push_surface_plot(&(plot_push_desc){
						.space = &sp,
						.data = z,
						.label = plot_snprintf("numerical %d state", n),
						.offset = res.eigenvalues[n],
						});
			}

			plot_update_until_closed();
			plot_shutdown();
#endif
		}

		sbmf_shutdown();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Solving hamiltonian for Harmonic oscillator basis functions

static i32 row_eigenfunction = 0;
static i32 col_eigenfunction = 0;

static f64 temp_potential(f64 x) {
	f64 edge = 3.0;
	f64 corner_radius = 0.5;
	if (fabs(x) > edge) {
		return 20.0;
	} else if (fabs(fabs(x) - edge) <= corner_radius) {
		f64 newx = fabs(x) - (edge-corner_radius);
		return corner_radius - sqrt(corner_radius*corner_radius - newx*newx);
	} else {
		return 0;
	}
}

static f64 compute_matrix_element(f64 x, void* data) {
	return (ho_eigenfunction(&row_eigenfunction, &x, 1)) * (temp_potential(x) - ho_potential(&x,1,0)) * (ho_eigenfunction(&col_eigenfunction, &x, 1));
}

describe(basis_ho) {
	it("1D"){
		sbmf_init();

		const i32 states_to_include = 16;
		const i32 matrix_order = states_to_include*states_to_include;

		c64 data[states_to_include*states_to_include];
		hermitian_bandmat mat = {
			.base = {
				.is_row_major = true,
				.rows = states_to_include,
				.cols = states_to_include,
				.data = data,
			},
			.bandcount = states_to_include,
			.size = states_to_include,
		};

		integration_settings settings = {
			.order = 7,
			.abs_error_tol = 1e-15,
			.rel_error_tol = 1e-1,
			.max_evals = 1e4,
		};

		integration_result res;

		printf("Matrix\n");
		for (i32 row = 0; row < states_to_include; ++row) {
			i32 starting_index = (states_to_include-1)*states_to_include + row;
			for (i32 col = row; col < states_to_include; ++col) {
				row_eigenfunction = row;
				col_eigenfunction = col;
				//																						(0,3) - 3
				//		(0,0)	(0,1)	(0,2)	(0,3)										(0,2) - 6	(1,3) - 7
				//		(1,0)	(1,1)	(1,2)	(1,3)		->					(0,1) - 9	(1,2) - 10	(2,3) - 11
				//		(2,0)	(2,1)	(2,2)	(2,3)				(0,0) - 12	(1,1) - 13	(2,2) - 14	(3,3) - 15
				//		(3,0)	(3,1)	(3,2)	(3,3)
				//0.246380 (y)    -0.000000 (y)   -0.290000 (y)   -0.000000 (y)   0.048792 (y)
				//-0.000000 (y)   0.851617 (y)    0.000000 (n)    -0.410967 (y)   0.000000 (n)
				//-0.290000 (y)   -0.000000 (n)   1.604115 (y)    0.000000 (y)    -0.521686 (y)
				//-0.000000 (y)   -0.410967 (y)   0.000000 (y)    2.244450 (y)    0.000000 (n)
				//0.048792 (y)    0.000000 (n)    -0.521686 (y)   -0.000000 (n)   3.360399 (y)

				//printf("(%d,%d) -- %d\t", row, col, starting_index - (states_to_include-1)*(col-row));
				res = quadgk(compute_matrix_element, -INFINITY,INFINITY, settings);
				////if (!res.converged)
				////	printf("matrix element %d,%d failed to converge!\n", row, col);
				////assert(res.converged);

				i32 out_index = starting_index - (states_to_include-1)*(col-row);
				c64* val = &data[out_index];
				*val = res.integral;

				if (row == col) {
					*val += ho_eigenvalue(&row, 1);
				}

				printf("%lf (%s)\t", creal(*val), res.converged ? "y" : "n");
			}
			printf("\n");
		}

		eig_result eig_res = eig_sparse_bandmat(mat, 3, EV_SMALLEST_RE);

		// Plotting
		{
#if PLOT_SPARSE_1D_PB_HOBASIS
			plot_init(800, 600, "1D PB HOBASIS");
			sample_space sp = make_linspace(1, -5,5, 150);
			grid g = generate_grid(1, (f64[]){-5}, (f64[]){5}, (i32[]){150});
			f32 z[sp.pointcount];

			FOREACH_POINT_1D(sp, i)
				z[i] = temp_potential(sp.points[i]);
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = z,
					.label = "potential",
					});

			for (u32 n = 0; n < eig_res.num_eigenpairs; ++n) {
				printf("Eigenvector: %d\n", n);
				memset(z, 0, sp.pointcount*sizeof(f32));
				for (u32 v = 0; v < eig_res.points_per_eigenvector; ++v) {
					printf("\t%lf\n", eig_res.eigenvectors[n*eig_res.points_per_eigenvector + v]);
					FOREACH_POINT_1D(sp, i)
						z[i] += eig_res.eigenvectors[n*eig_res.points_per_eigenvector +  v]*ho_eigenfunction(&v, (f64*)&g.points[i], 1);
				}
				push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = z,
					.label = plot_snprintf("numerical %d state", n),
					.offset = eig_res.eigenvalues[n],
					});
			}

			plot_update_until_closed();
			plot_shutdown();
#endif
		}

		sbmf_shutdown();
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Basic operations

describe(matrix_ops) {
	// Produced by matlab, (warning COLUMN MAJOR!)
	static c64 eigvecs_answer[] = {
		-0.4489,    0.0085,    0.2281,   -0.2686,    0.3344,    0.4139,    0.2168,   -0.4228,   -0.3533,    0.2015,
		0.6377 ,  -0.0089 ,  -0.2163 ,   0.1788 ,  -0.1397 ,   0.0469 ,   0.0917 ,  -0.3495 ,  -0.4866 ,   0.3542 ,
		-0.5494,   -0.0178,   -0.0802,    0.2295,   -0.3348,   -0.4370,   -0.1732,    0.1008,   -0.3448,    0.4203,
		0.2818 ,   0.1498 ,   0.4358 ,  -0.4890 ,   0.1939 ,  -0.1808 ,  -0.2166 ,   0.3980 ,  -0.0938 ,   0.4297 ,
		-0.0251,   -0.4047,   -0.3178,    0.3091,    0.5092,    0.2221,    0.0615,    0.3800,    0.1375,    0.4042,
		-0.0817,    0.3838,   -0.3454,   -0.2371,   -0.4336,    0.4985,   -0.0242,    0.1086,    0.3041,    0.3595,
		0.0342 ,   0.2572 ,   0.4964 ,   0.4522 ,  -0.0600 ,  -0.0767 ,   0.5003 ,  -0.0762 ,   0.3593 ,   0.2951 ,
		0.0375 ,  -0.6138 ,   0.2489 ,  -0.0935 ,  -0.2515 ,   0.0495 ,  -0.3502 ,  -0.3963 ,   0.3852 ,   0.2395 ,
		-0.0204,   -0.0829,   -0.3886,   -0.4527,    0.1130,   -0.5182,    0.4676,   -0.2221,    0.2552,    0.1495,
		-0.0189,    0.4651,   -0.1586,    0.1849,    0.4382,   -0.1679,   -0.5212,   -0.4027,    0.2381,    0.1154
	};

	static f64 eigvals_answer[] = {
		-11.785562861423848,
		-8.41709283554329,
		-7.5370741723957915,
		-4.989284977718224,
		-2.7593112024185578,
		2.0207266054245863,
		4.806575870197963,
		8.440803849012717,
		13.3962999372205,
		16.823919787643952,
	};

	before_each() {
		sbmf_init();
	}
	after_each() {
		sbmf_shutdown();
	}

	it ("dense (symetric) upper tridiagonal band matrix") {
		c64 bands[] = {
			0,0,0,1,2,3,4,5,6,7,
			0,9,8,7,6,5,4,3,2,1,
			1,1,1,1,1,1,1,1,1,1,
		};

		i32 size = 10;

		hermitian_bandmat bm = {
			.base = {
				.is_row_major = true,
				.rows = 3,
				.cols = size,
				.data = bands
			},
			.bandcount = 3,
			.size = size,
		};

		f64 eigvals[size];
		c64 eigvecs[size*size];
		// Only use upper portion (including main diag) of bandmat bm.
		// Also, this solves for all eigenvalues!
		eig_dense_symetric_upper_tridiag_bandmat(bm, eigvals, eigvecs);
		// check eigenvalues
		for (i32 i = 0; i < size; ++i) {
			assert(fabs(eigvals[i] - eigvals_answer[i]) <= 0.01);
		}

		// Check eigenvectors
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				c64 got = eigvecs[i*size + j];
				c64 expected = eigvecs_answer[j*size + i];
				assert(cabs(got)-cabs(expected) <= 1);
			}
		}
	}

	it ("sparse bandmatrix") {
		c64 bands[] = {
			0,0,0,1,2,3,4,5,6,7,
			0,9,8,7,6,5,4,3,2,1,
			1,1,1,1,1,1,1,1,1,1,
		};

		i32 size = 10;

		hermitian_bandmat bm = {
			.base = {
				.is_row_major = true,
				.rows = 3,
				.cols = size,
				.data = bands
			},
			.bandcount = 3,
			.size = size,
		};

		u32 eigenpairs_to_find = 4;
		eig_result res = eig_sparse_bandmat(bm, eigenpairs_to_find, EV_SMALLEST_RE);

		// check eigenvalues
		for (u32 i = 0; i < res.num_eigenpairs; ++i) {
			assert(cabs(res.eigenvalues[i] - eigvals_answer[i]) <= 0.001);
		}

		// Check eigenvectors
		for (u32 i = 0; i < res.num_eigenpairs; ++i) {
			for (u32 j = 0; j < res.points_per_eigenvector; ++j) {
				c64 got = res.eigenvectors[i*bm.size + j];
				c64 expected = eigvecs_answer[j*bm.size + i];
				assert(cabs(got)-cabs(expected) <= 0.01);
			}
		}
	}

	it ("hermitian bandmat mulv") {
		c64 data[10] = {
			0,0,1,0,0,
			1,1,1,1,1,
		};

		hermitian_bandmat m = {
			.base = {
				.is_row_major = true,
				.rows = 2,
				.cols = 5,
				.data = data,
			},
			.bandcount = 2,
			.size = 5,
		};

		c64 input[5] = {1,1,1,1,1};
		c64 ans[5];
		c64 expected[5] = {1,2,2,1,1};
		complex_hermitian_bandmat_mulv(ans, m, input);

		for (i32 i = 0; i < 5; ++i) {
			asserteq(ans[i], expected[i]);
		}
	}
}

describe(simple_math_funcs) {
	it("hermite polynomials") {
		asserteq(hermite_poly(0,0), 1);
		asserteq(hermite_poly(1,5), 2*5);
		asserteq(hermite_poly(2,5), 4*5*5-2);
	}

	it("factorial") {
		asserteq(factorial(0), 1);
		asserteq(factorial(5), 120);
		asserteq(factorial(10), 3628800);
		asserteq(factorial(20), 2432902008176640000);
	}
}

snow_main();
