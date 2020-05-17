#define SNOW_ENABLED
#include "snow.h"

#include <sbmf/groundstate_solver/groundstate_solver.h>
#include <sbmf/common/profile.h>
#include <sbmf/common/fdm.h>
#include <sbmf/common/grid.h>

#include <sbmf/common/eigenproblem.h>

//#include <stdio.h>
//#include <unistd.h>

#include <plot/plot.h>

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
// TEST SOLVING HAMILTONIAN VIA HARMONIC OSC. BASIS FUNCTIONS
#include <sbmf/basis/harmonic_oscillator.h>

static i32 row_eigenfunction = 0;
static i32 col_eigenfunction = 0;

static f64 temp_potential(f64 x) {
	return ho_potential(&x, 1, 0);
}

static f64 compute_matrix_element(f64 x, void* data) {
	return (ho_eigenfunction(&row_eigenfunction, &x, 1)) * (temp_potential(x) - ho_potential(&x, 1, 0)) * (ho_eigenfunction(&col_eigenfunction, &x, 1));
}

describe(solve_hamiltonian_ho_basis) {
	it("temp"){
		static const i32 states_to_include = 5;
		static const i32 matrix_order = states_to_include*states_to_include;

		f64 mat[matrix_order];
		memset(mat, 0, matrix_order*sizeof(i32));

		integration_settings settings = {
			.order = 7,
			.abs_error_tol = 0,
			.rel_error_tol = 1e-1,
			.max_evals = 1e4,
		};

		integration_result res;

		for (i32 row = 0; row < states_to_include; ++row) {
			for (i32 col = 0; col < states_to_include; ++col) {
				row_eigenfunction = row;
				col_eigenfunction = col;

				res = quadgk(compute_matrix_element, -INFINITY,INFINITY, settings);
				assert(res.converged);

				mat[row*states_to_include + col] = res.integral;

				if (row == col) {
					mat[col*states_to_include + col] += ho_eigenvalue(&row, 1);
				}
			}
		}

	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// TEST GRID BASED EIGENVALUE SOLVER AGAINST KNOWN SOLUTIONS

static void normalize_function_c64(c64* func, i32 size) {
	f64 sum = 0.0;
	for (i32 i = 0; i < size; ++i) {
		sum += cabs(func[i])*cabs(func[i]);
	}

	f64 scaling_factor = 1.0/sqrt(sum);
	for (i32 i = 0; i < size; ++i) {
		func[i] *= scaling_factor;
	}
}
static void normalize_function_f64(f64* func, i32 size) {
	f64 sum = 0.0;
	for (i32 i = 0; i < size; ++i) {
		sum += func[i]*func[i];
	}

	f64 scaling_factor = 1.0/sqrt(sum);
	for (i32 i = 0; i < size; ++i) {
		func[i] *= scaling_factor;
	}
}
static void normalize_function_f32(f32* func, i32 size) {
	f32 sum = 0.0;
	for (i32 i = 0; i < size; ++i) {
		sum += func[i]*func[i];
	}

	f32 scaling_factor = 1.0/sqrt(sum);
	for (i32 i = 0; i < size; ++i) {
		func[i] *= scaling_factor;
	}
}

#define normalize_function(f, size) 									\
	_Generic((f), 																			\
			c64*: normalize_function_c64(f, size),					\
			f64*: normalize_function_f64(f, size),					\
			f32*: normalize_function_f32(f, size),					\
			default: assert(0))

f64 periodic_pot(f64* v, i32 n, c64 u) {
	f64 temp = 0.0;
	for (i32 i = 0; i < n; ++i) {
		temp += cos(v[i])*cos(v[i]);
	}
	return temp;
}

c64 initial_guess(f64* v, i32 n) {
	return 1.0/10*10;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
describe(grid_based_ho_hamiltonian_solving) {
	return;
	const f64 assert_limit = 0.1;

	it ("1D") {
		const i32 N = 64;

		grid g = generate_grid(1, 
				(f64[]){-10},
				(f64[]){ 10},
				(i32[]){ N}
				);

		bandmat fdm;
		{
			// NOTE: g.pointcounts[0] forces it to be square!
			fdm = generate_fd_matrix(g.pointcounts[0], g.dimensions, g.deltas);

			for (i32 i = 0; i < fdm.size*(fdm.super_diags+1); ++i) {
				fdm.bands[i] = -0.5*fdm.bands[i];
			}

			for (i32 i = 0; i < fdm.size; ++i) {
				i32 idx = fdm.size*(fdm.super_diags) + i;
				fdm.bands[idx] += (c64) ho_potential(&g.points[i], g.dimensions, 0.0);
			}	
		}

		f64* eigenvalues = malloc(fdm.size * sizeof(f64));
		c64* eigenvectors = malloc((fdm.size*fdm.size) * sizeof(c64));
		eig_dense_symetric_upper_tridiag_bandmat(fdm, eigenvalues, eigenvectors);

		f64 expected_answer[g.total_pointcount];
		for (i32 i = 0; i < 10; ++i) {
			normalize_function_c64(&eigenvectors[i*g.total_pointcount], g.total_pointcount);

			for (i32 j = 0; j < g.total_pointcount; ++j) {
				expected_answer[j] = ho_eigenfunction(&i, &g.points[j], g.dimensions);
			}
			normalize_function_f64(expected_answer, g.total_pointcount);

			for (i32 j = 0; j < g.total_pointcount; ++j) {
				f64 expected = fabs(expected_answer[j]);
				f64 got = fabs(creal(eigenvectors[i*g.total_pointcount + j]));
				assert(expected-got < assert_limit);
			}

			f64 expected_energy = ho_eigenvalue(&i, g.dimensions);
			assert(eigenvalues[i] - expected_energy < assert_limit);
		}

		(void)normalize_function_f32;
		free(fdm.bands);
		free(eigenvalues);
		free(eigenvectors);
		free_grid(g);
	}
	it ("2D") {
		const i32 N = 64;

		grid g = generate_grid(2, 
				(f64[]){-10,-10},
				(f64[]){ 10, 10},
				(i32[]){ N,  N}
				);

		bandmat fdm;
		{
			// NOTE: g.pointcounts[0] forces it to be square!
			fdm = generate_fd_matrix(g.pointcounts[0], g.dimensions, g.deltas);

			for (i32 i = 0; i < fdm.size*(fdm.super_diags+1); ++i) {
				fdm.bands[i] = -0.5*fdm.bands[i];
			}

			for (i32 i = 0; i < fdm.size; ++i) {
				i32 idx = fdm.size*(fdm.super_diags) + i;
				fdm.bands[idx] += (c64) ho_potential(&g.points[i], g.dimensions, 0.0);
			}	
		}

		f64* eigenvalues = malloc(fdm.size * sizeof(f64));
		c64* eigenvectors = malloc((fdm.size*fdm.size) * sizeof(c64));
		eig_dense_symetric_upper_tridiag_bandmat(fdm, eigenvalues, eigenvectors);

		f64 expected_answer[g.total_pointcount];
		for (i32 i = 0; i < 10; ++i) {
			normalize_function_c64(&eigenvectors[i*g.total_pointcount], g.total_pointcount);

			for (i32 j = 0; j < g.total_pointcount; ++j) {
				expected_answer[j] = ho_eigenfunction(&i, &g.points[j], g.dimensions);
			}
			normalize_function_f64(expected_answer, g.total_pointcount);

			for (i32 j = 0; j < g.total_pointcount; ++j) {
				f64 expected = fabs(expected_answer[j]);
				f64 got = fabs(creal(eigenvectors[i*g.total_pointcount + j]));
				assert(expected-got < assert_limit);
			}

			f64 expected_energy = ho_eigenvalue(&i, g.dimensions);
			assert(eigenvalues[i] - expected_energy < assert_limit);
		}

		(void)normalize_function_f32;
		free(fdm.bands);
		free(eigenvalues);
		free(eigenvectors);
		free_grid(g);
	}
	it ("3D") {
		const i32 N = 64;

		grid g = generate_grid(3, 
				(f64[]){-10,-10,-10},
				(f64[]){ 10, 10, 10},
				(i32[]){ N,  N,  N}
				);

		bandmat fdm;
		{
			// NOTE: g.pointcounts[0] forces it to be square!
			fdm = generate_fd_matrix(g.pointcounts[0], g.dimensions, g.deltas);

			for (i32 i = 0; i < fdm.size*(fdm.super_diags+1); ++i) {
				fdm.bands[i] = -0.5*fdm.bands[i];
			}

			for (i32 i = 0; i < fdm.size; ++i) {
				i32 idx = fdm.size*(fdm.super_diags) + i;
				fdm.bands[idx] += (c64) ho_potential(&g.points[i], g.dimensions, 0.0);
			}	
		}

		f64* eigenvalues = malloc(fdm.size * sizeof(f64));
		c64* eigenvectors = malloc((fdm.size*fdm.size) * sizeof(c64));
		eig_dense_symetric_upper_tridiag_bandmat(fdm, eigenvalues, eigenvectors);

		f64 expected_answer[g.total_pointcount];
		for (i32 i = 0; i < 10; ++i) {
			normalize_function_c64(&eigenvectors[i*g.total_pointcount], g.total_pointcount);

			for (i32 j = 0; j < g.total_pointcount; ++j) {
				expected_answer[j] = ho_eigenfunction(&i, &g.points[j], g.dimensions);
			}
			normalize_function_f64(expected_answer, g.total_pointcount);

			for (i32 j = 0; j < g.total_pointcount; ++j) {
				f64 expected = fabs(expected_answer[j]);
				f64 got = fabs(creal(eigenvectors[i*g.total_pointcount + j]));
				assert(expected-got < assert_limit);
			}

			f64 expected_energy = ho_eigenvalue(&i, g.dimensions);
			assert(eigenvalues[i] - expected_energy < assert_limit);
		}

		(void)normalize_function_f32;
		free(fdm.bands);
		free(eigenvalues);
		free(eigenvectors);
		free_grid(g);
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline void test_eigenvalue_solving_harmonic_osc_3d() { }
static inline void test_eigenvalue_solving_particle_in_a_box_1d() { }
static inline void test_eigenvalue_solving_particle_in_a_box_2d() { }
static inline void test_eigenvalue_solving_particle_in_a_box_3d() { }
static inline void test_eigenvalue_solving_periodic_pot_1d() {
	//printf("-- Running 1D periodic test\n");

	const i32 N = 256;

	grid g = generate_grid(1, 
			(f64[]){-5},
			(f64[]){ 5},
			 (i32[]){ N}
			);

	PROFILE_BEGIN("gen. fdm ho 1d");
	bandmat fdm;
	{
		// NOTE: g.pointcounts[0] forces it to be square!
		fdm = generate_fd_matrix(g.pointcounts[0], g.dimensions, g.deltas);

		for (i32 i = 0; i < fdm.size*(fdm.super_diags+1); ++i) {
			fdm.bands[i] = -0.5*fdm.bands[i];
		}

		for (i32 i = 0; i < fdm.size; ++i) {
			i32 idx = fdm.size*(fdm.super_diags) + i;
			fdm.bands[idx] += (c64) periodic_pot(&g.points[i], g.dimensions, 0.0);
		}	
	}
	PROFILE_END("gen. fdm  ho 1d");


	PROFILE_BEGIN("e.v. prob. ho 1d");
	f64* eigenvalues = malloc(fdm.size * sizeof(f64));
	c64* eigenvectors = malloc((fdm.size*fdm.size) * sizeof(c64));
	eig_dense_symetric_upper_tridiag_bandmat(fdm, eigenvalues, eigenvectors);
	PROFILE_END("e.v. prob. ho 1d");


	plotstate* state = plot_init();

	f32 x[g.total_pointcount];
	f32 v[g.total_pointcount];
	f32 u[g.total_pointcount];
	for (i32 i = 0; i < g.total_pointcount; ++i) {
		x[i] = g.points[i];
		v[i] = periodic_pot(&g.points[i], g.dimensions, 0.0);
	}

	plot_1d(state, x, v, g.total_pointcount);

	for (i32 i = 0; i < 5; ++i) {
		normalize_function_c64(&eigenvectors[i*g.total_pointcount], g.total_pointcount);
		for (i32 j = 0; j < g.total_pointcount; ++j) {
			u[j] = 10*creal(eigenvectors[i*g.total_pointcount + j]);
		}

		plot_1d(state, x, u, g.total_pointcount);
	}


	plot_wait_on_join(state);
	plot_shutdown(state);

	(void)normalize_function_f32;
	free_grid(g);
}

describe(solve_eigenproblems_bandmat) {
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

	it ("dense (symetric) upper tridiagonal band matrix") {
		c64 bands[] = {
			0,0,0,1,2,3,4,5,6,7,
			0,9,8,7,6,5,4,3,2,1,
			1,1,1,1,1,1,1,1,1,1,
			9,8,7,6,5,4,3,2,1,0,
			0,1,2,3,4,5,6,7,0,0,
		};

		bandmat bm = {
			.size = 10,
			.super_diags = 2,
			.sub_diags = 2,
			.bands = bands,
		};

		f64 eigvals[bm.size];
		c64 eigvecs[bm.size*bm.size];
		// Only use upper portion (including main diag) of bandmat bm.
		// Also, this solves for all eigenvalues!
		eig_dense_symetric_upper_tridiag_bandmat(bm, eigvals, eigvecs);
		// check eigenvalues
		for (i32 i = 0; i < bm.size; ++i) {
			assert(fabs(eigvals[i] - eigvals_answer[i]) <= 1);
		}

		// Check eigenvectors
		for (int i = 0; i < bm.size; ++i) {
			for (int j = 0; j < bm.size; ++j) {
				c64 got = eigvecs[i*bm.size + j]; 
				c64 expected = eigvecs_answer[j*bm.size + i];
				assert(cabs(got)-cabs(expected) <= 1);
			}
		}
	}

	it ("sparse bandmatrix") {
		c64 bands[] = {
			0,0,0,1,2,3,4,5,6,7,
			0,9,8,7,6,5,4,3,2,1,
			1,1,1,1,1,1,1,1,1,1,
			9,8,7,6,5,4,3,2,1,0,
			0,1,2,3,4,5,6,7,0,0,
		};

		bandmat bm = {
			.size = 10,
			.super_diags = 2,
			.sub_diags = 2,
			.bands = bands,
		};

		c64 eigvals[bm.size];
		c64 eigvecs[bm.size*bm.size];

		i32 eigenpairs_to_find = 4;
		eig_sparse_bandmat(bm, eigenpairs_to_find, EV_SMALLEST_RE, eigvals, eigvecs);

		//// check eigenvalues
		//for (i32 i = 0; i < eigenpairs_to_find; ++i) {
		//	printf("%lf + i%lf\n", creal(eigvals[i]), cimag(eigvals[i]));
		//	//assert(fabs(eigvals[i] - eigvals_answer[i]) <= 1);
		//}

		//// Check eigenvectors
		//for (int i = 0; i < eigenpairs_to_find; ++i) {
		//	for (int j = 0; j < bm.size; ++j) {
		//		c64 got = eigvecs[i*bm.size + j]; 
		//		c64 expected = eigvecs_answer[j*bm.size + i];
		//		//assert(cabs(got)-cabs(expected) <= 1);
		//	}
		//}
	}

	it ("trivial bandmatrix example 1") {
		bandmat bm = {
			.size = 3,
			.super_diags = 1,
			.sub_diags = 1,
			.bands = (c64*)(c64[]){
				0,4,5,
				1,2,3,
				4,5,0,
			},
		};

		c64 eigvals[bm.size];
		c64 eigvecs[bm.size*bm.size];

		i32 eigenpairs_to_find = 1;
		eig_sparse_bandmat(bm, eigenpairs_to_find, EV_SMALLEST_MAG, eigvals, eigvecs);
	}

	it ("trivial bandmatrix example 2") {
		bandmat bm = {
			.size = 3,
			.super_diags = 2,
			.sub_diags = 2,
			.bands = (c64*)(c64[]){
				0,0,1,
				0,0,0,
				1,2,3,
				0,0,0,
				1,0,0
			},
		};

		c64 eigvals[bm.size];
		c64 eigvecs[bm.size*bm.size];

		i32 eigenpairs_to_find = 1;
		eig_sparse_bandmat(bm, eigenpairs_to_find, EV_SMALLEST_MAG, eigvals, eigvecs);
	}


	it ("trivial real matrix example") {
		f64 mat[3*3] = {
			1,4,0,
			4,2,5,
			0,5,3
		};

		eig_sparse_real(mat, 3, 1, NULL,NULL);
	}
}

#include <cblas.h>
describe(test_cblas) {
	it ("zgbmv") {
		i32 kl = 2, ku = 2;
		i32 n = 10;
		i32 lda = 5;

		c64 bands[] = {
			0,0,0,1,2,3,4,5,6,7,
			0,9,8,7,6,5,4,3,2,1,
			1,1,1,1,1,1,1,1,1,1,
			9,8,7,6,5,4,3,2,1,0,
			0,1,2,3,4,5,6,7,0,0,
		};

		c64 input_mat[lda*n];
		mat_transpose(input_mat, bands, lda, n);

		c64 x[10] = {0,1,2,3,4,5,6,7,8,9};
		c64 y[n];

		c64 zero = 0, one = 1;
		cblas_zgbmv(CblasColMajor, CblasNoTrans, n, n, kl, ku, &one, input_mat, 
				lda, x, 1, &zero, y, 1);

		//	1 9 0                   0
		//	9 1 8 1 								1
		// 	0 8 1 7 2 							2
		//    1 7 1 6 3 						3
		//      2 6 1 5 4 					4
		//        3 5 1 4 5 				5
		//          4 4 1 3 6 			6
		//            5 3 1 2 7			7
		//              6 2 1 1			8
		//                7 1 1			9
		
		//	9                   
		//	20 								
		// 	39
		//  57
		//  72
		//  93
		//  111
		//  129
		//  67
		//  66
		
		c64 expected_y[10] = {9,20,39,57,75,93,111,129,67,66};
		for (i32 i = 0; i < n; ++i) {
			asserteq(expected_y[i], y[i]);
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

describe(common) {
	it ("matrix transpose") {
		i32 m = 5;
		i32 n = 3;
		c64 mat[5*3] = {
			1,  2, 3,
			4,  5, 6,
			8,  9,10,
			11,12,13,
			14,15,16
		};

		c64 ans[m*n];
		mat_transpose(ans, mat, m,n);
		mat_transpose(ans, ans, n,m);

		for (i32 i = 0; i < m*n; ++i) {
			asserteq(ans[i], mat[i]);
		}
	}
}

snow_main();
