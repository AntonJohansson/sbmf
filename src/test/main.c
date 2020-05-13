#define SNOW_ENABLED
#include "snow.h"

#include <sbmf/groundstate_solver/groundstate_solver.h>
#include <sbmf/common/profile.h>
#include <sbmf/common/fdm.h>
#include <sbmf/common/grid.h>
#include <sbmf/common/eigenproblem.h>

#include <sbmf/basis/harmonic_oscillator.h>

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

// n-dimensional harmonic osc potential.
f64 harmonic_osc_potential(f64* v, i32 n, c64 u) {
	f64 temp = 0.0;
	for (i32 i = 0; i < n; ++i) {
		temp += v[i]*v[i];
	}
	return temp*0.5;
}

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

static inline void test_eigenvalue_solving_harmonic_osc_1d() {
	const i32 N = 256;

	grid g = generate_grid(1, 
			(f64[]){-10},
			(f64[]){ 10},
			 (i32[]){ N}
			);

	PROFILE_BEGIN("gen. fdm ho 1d");
	bandmat fdm;
	{
		// NOTE: g.pointcounts[0] forces it to be square!
		fdm = generate_fd_matrix(g.pointcounts[0], g.dimensions, g.deltas);

		for (i32 i = 0; i < fdm.size*fdm.bandcount; ++i) {
			fdm.bands[i] = -0.5*fdm.bands[i];
		}

		for (i32 i = 0; i < fdm.size; ++i) {
			i32 idx = fdm.size*(fdm.bandcount-1) + i;
			fdm.bands[idx] += (c64) ho_potential(&g.points[i], g.dimensions, 0.0);
		}	
	}
	PROFILE_END("gen. fdm  ho 1d");


	PROFILE_BEGIN("e.v. prob. ho 1d");
	f64* eigenvalues = malloc(fdm.size * sizeof(f64));
	c64* eigenvectors = malloc(fdm.order * sizeof(c64));
	eigvp_bandmat(eigenvalues, eigenvectors, fdm);
	PROFILE_END("e.v. prob. ho 1d");

	plotstate* state = plot_init();
	f32 x[g.total_pointcount];
	f32 v[g.total_pointcount];
	for (i32 i = 0; i < g.total_pointcount; ++i) {
		x[i] = g.points[i];
		v[i] = ho_potential(&g.points[i], g.dimensions, 0.0);
	}
	plot_1d(state, x, v, g.total_pointcount);

	f64 expected_answer[g.total_pointcount];
	for (i32 i = 0; i < 10; ++i) {
		normalize_function_c64(&eigenvectors[i*g.total_pointcount], g.total_pointcount);
		
		for (i32 j = 0; j < g.total_pointcount; ++j) {
			//expected_answer[j] = psi_ho_1d(i, g.points[j]);
			expected_answer[j] = ho_eigenfunction(&i, &g.points[j], g.dimensions);
		}
		normalize_function_f64(expected_answer, g.total_pointcount);

		for (i32 j = 0; j < g.total_pointcount; ++j) {
			f64 expected = fabs(expected_answer[j]);
			f64 got = fabs(creal(eigenvectors[i*g.total_pointcount + j]));
			if (expected-got > 0.05) {
				//fprintf(stderr, "\t- Eigenvector calculation failed: got %lf; expected: %lf\n", got, expected);
			}
		}

		f64 expected_energy = ho_eigenvalue(&i, g.dimensions);
		if (fabs(eigenvalues[i] - expected_energy) > 0.05) {
			//fprintf(stderr, "\t - Eigenvalue calculation failed: got %lf, expected %lf\n", eigenvalues[i], expected_energy);
		}

		f32 y[g.total_pointcount];
		f32 z[g.total_pointcount];
		for (i32 j = 0; j < g.total_pointcount; ++j) {
			//y[j] = E_ho_1d(i) + expected_answer[j];
			y[j] = ho_eigenvalue(&i, g.dimensions) + expected_answer[j];
			z[j] = eigenvalues[i] + creal(eigenvectors[i*g.total_pointcount + j]);
		}

		plot_1d(state, x, y, g.total_pointcount);
		plot_1d(state, x, z, g.total_pointcount);
	}

	plot_wait_on_join(state);
	plot_shutdown(state);

	(void)normalize_function_f32;
	free_grid(g);
}

static inline void test_eigenvalue_solving_harmonic_osc_2d() {
	//printf("-- Running 2D harmonic osc test\n");

	const i32 N = 32;

	grid g = generate_grid(2, 
			(f64[]){-2, -2},
			(f64[]){ 2,  2},
			 (i32[]){ N,   N}
			);








	PROFILE_BEGIN("gen. fdm ho 2d");
	bandmat fdm;
	{
		// NOTE: g.pointcounts[0] forces it to be square!
		fdm = generate_fd_matrix(g.pointcounts[0], g.dimensions, g.deltas);

		for (i32 i = 0; i < fdm.size*fdm.bandcount; ++i) {
			fdm.bands[i] = -0.5*fdm.bands[i];
		}

		for (i32 i = 0; i < fdm.size; ++i) {
			i32 idx = fdm.size*(fdm.bandcount-1) + i;
			fdm.bands[idx] += (c64) harmonic_osc_potential(&g.points[2*i], g.dimensions, 0.0);
		}	
	}
	PROFILE_END("gen. fdm ho 2d");


	PROFILE_BEGIN("e.v. prob. ho 2d");
	f64* eigenvalues = malloc(fdm.size * sizeof(f64));
	c64* eigenvectors = malloc(fdm.order * sizeof(c64));
	eigvp_bandmat(eigenvalues, eigenvectors, fdm);
	PROFILE_END("e.v. prob. ho 2d");

	normalize_function_c64(&eigenvectors[0], g.total_pointcount);

	f32 x[g.total_pointcount];
	f32 y[g.total_pointcount];
	f32 v[g.total_pointcount];
	f32 u[g.total_pointcount];
	f32 e[g.total_pointcount];
	for (i32 i = 0; i < g.total_pointcount; ++i) {
		x[i] = g.points[2*i];
		y[i] = g.points[2*i+1];
		v[i] = ho_potential(&g.points[2*i], g.dimensions, 0.0);
		u[i] = -10*creal(eigenvectors[1*g.total_pointcount + i]);
		e[i] = ho_eigenfunction((i32[]){1,0}, &g.points[2*i], g.dimensions);
	}

	for (i32 i = 0; i < g.total_pointcount; ++i) {
		if (fabs(e[i] - u[i]) > 0.01) {
			//fprintf(stderr, "\t- Eigenvector calculation %d failed: got %lf; expected %lf\n", i, e[i], u[i]);
		}
	}

	normalize_function_f32(&e[0], g.total_pointcount);
	for (i32 i = 0; i < g.total_pointcount; ++i) {
		e[i] = 10*e[i];
	}

	plotstate* state = plot_init();
	plot_2d(state, x, y, v, g.total_pointcount);
	plot_2d(state, x, y, u, g.total_pointcount);
	plot_2d(state, x, y, e, g.total_pointcount);
	plot_wait_on_join(state);
	plot_shutdown(state);
	
	
	
	
	

	
	
	
	

	free(eigenvalues);
	free(eigenvectors);
	free_grid(g);
}

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

		for (i32 i = 0; i < fdm.size*fdm.bandcount; ++i) {
			fdm.bands[i] = -0.5*fdm.bands[i];
		}

		for (i32 i = 0; i < fdm.size; ++i) {
			i32 idx = fdm.size*(fdm.bandcount-1) + i;
			fdm.bands[idx] += (c64) periodic_pot(&g.points[i], g.dimensions, 0.0);
		}	
	}
	PROFILE_END("gen. fdm  ho 1d");


	PROFILE_BEGIN("e.v. prob. ho 1d");
	f64* eigenvalues = malloc(fdm.size * sizeof(f64));
	c64* eigenvectors = malloc(fdm.order * sizeof(c64));
	eigvp_bandmat(eigenvalues, eigenvectors, fdm);
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

////////////////////////////////////////////////////////////////////////////////////////////////////
// TEST GRID BASED EIGENVALUE SOLVER AGAINST KNOWN SOLUTIONS

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

describe(bandmat_eigenproblem) {
	it("basic solving") {
		c64 bands[] = {
			0,0,0,1,2,3,4,5,6,7,
			0,9,8,7,6,5,4,3,2,1,
			1,1,1,1,1,1,1,1,1,1
		};

		bandmat bm = {
			.order = 100,
			.size = 10,
			.bandcount = 3,
			.bands = bands,
		};

		f64 eigvals[bm.size];
		c64 eigvecs[bm.order];
		eigvp_bandmat(eigvals, eigvecs, bm);

		// Produced by matlab, (warning COLUMN MAJOR!)
		c64 eigvecs_answer[] = {
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

		f64 eigvals_answer[] = {
			-11.7856,
			-8.4171 ,
			-7.5371 ,
			-4.9893 ,
			-2.7593 ,
			2.27    ,
			4.866   ,
			8.448   ,
			13.3963 ,
			16.8239 
		};

		// check eigenvalues
		for (i32 i = 0; i < bm.size; ++i) {
			assert(fabs(eigvals[i] - eigvals_answer[i]) <= 1);
			//if (fabs(eigvals[i] - eigvals_answer[i]) > 0.01) {
			//	fprintf(stderr, "\t- Eigenvalue calculation %d failed: got %lf; expected %lf\n", i, eigvals[i], eigvals_answer[i]);
			//}
		}

		// Check eigenvectors
		for (int i = 0; i < bm.size; ++i) {
			for (int j = 0; j < bm.size; ++j) {
				c64 got = eigvecs[i*bm.size + j]; 
				c64 expected = eigvecs_answer[j*bm.size + i];
				assert(cabs(got)-cabs(expected) <= 1);
				//if (cabs(got)-cabs(expected) > 0.01) {
				//	fprintf(stderr, "\t- Eigenvector calculation failed: got %lf + %lf i; expected: %lf + %lf i\n", creal(got), cimag(got), creal(expected), cimag(expected));
				//}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// TEST GRID BASED EIGENVALUE SOLVER AGAINST KNOWN SOLUTIONS

//int main(int argc, char** argv) {
//	assert(hermite_poly(0,0) == 1);
//	assert(hermite_poly(1,5) == 2*5);
//	assert(hermite_poly(2,5) == 4*5*5-2);
//
//	assert(factorial(0) == 1);
//	assert(factorial(5) == 120);
//	assert(factorial(10) == 3628800);
//	assert(factorial(20) == 2432902008176640000);
//
//	test_prioqueue();
//	test_quadgk();
//
//	return 0;
//
//	printf("Testing eigenvalue solving\n"),
//	//test_eigenvalue_solving();
//	test_eigenvalue_solving_harmonic_osc_1d();
//	test_eigenvalue_solving_periodic_pot_1d();
//	//test_eigenvalue_solving_harmonic_osc_2d();
//
//	//plotstate* state = plot_init();
//
//	//const i32 N = 128;
//
//	//grid g = generate_grid(1, 
//	//		(f64[]){-5, -5},
//	//		(f64[]){ 5,  5},
//	//		 (i32[]){ N,  	 	 N     }
//	//		);
//
//	//gss_settings settings = {
//	//	.g = g,
//	//	.max_iterations = 2000,
//	//	.error_tol = 1e-10,
//	//	.dt = 0.01,
//
//	//	.measure_every = 0,
//	//};
//
//	//PROFILE_BEGIN("item");
//	//gss_result res = item_execute(settings,
//	//															potential,
//	//															initial_guess);
//	//PROFILE_END("item");
//	//printf(
//	//		"ITEM results:\n"
//	//		"\titerations: %d/%d\n"
//	//		"\terror: %e\n"
//	//		"\terror tol: %e\n",
//	//	res.iterations,
//	//	res.settings.max_iterations,
//	//	res.error,
//	//	res.settings.error_tol
//	//);
//
//
//	//float x[g.total_pointcount];
//	//float y[g.total_pointcount];
//	//float v[g.total_pointcount];
//	//for (i32 i = 0; i < g.total_pointcount; ++i) {
//	//	x[i] = (float)g.points[i];
//	//	y[i] = (float)creal(res.wavefunction[i]);
//	//	v[i] = (float)potential(&g.points[i], 1, res.wavefunction[i]);
//	//}
//
//	//plot_1d(state, x, v, g.total_pointcount); 
//	//plot_1d(state, x, y, g.total_pointcount); 
//
//
//	//PROFILE_BEGIN("gen. fdm");
//	//bandmat fdm;
//	//{
//	//	// NOTE: g.pointcounts[0] forces it to be square!
//	//	fdm = generate_fd_matrix(g.pointcounts[0], pow(2,g.dimensions), g.dimensions, g.deltas);
//
//	//	for (i32 i = 0; i < fdm.size*fdm.bandcount; ++i) {
//	//		fdm.bands[i] = -0.5*fdm.bands[i];
//	//	}
//
//	//	for (i32 i = 0; i < fdm.size; ++i) {
//	//		i32 idx = fdm.size*(fdm.bandcount-1) + i;
//	//		fdm.bands[idx] += (c64) potential(&g.points[i], 1, res.wavefunction[i]);
//	//	}	
//	//}
//	//PROFILE_END("gen. fdm");
//
//
//	//PROFILE_BEGIN("e.v. prob.");
//	//{
//	//	f64* eigenvalues = malloc(fdm.size * sizeof(f64));
//	//	c64* eigenvectors = malloc(fdm.order * sizeof(c64));
//
//	//	eigvp_bandmat(eigenvalues, eigenvectors, fdm);
//
//	//	float x[fdm.size];
//	//	float u[fdm.size];
//	//	for (int i = 1; i < 6; ++i) {
//	//		for (int j = 0; j < fdm.size; ++j) {
//	//			x[j] = (float)g.points[j];
//	//			u[j] = (float)eigenvalues[i] +  3*(float)creal(eigenvectors[j*fdm.size + i]);
//	//		}
//	//		plot_1d(state, x,u, fdm.size);
//	//	}
//	//}
//	//PROFILE_END("e.v. prob.");
//
//
//
//
//	//gss_free_result(res);
//	//free_grid(g);
//
//	profile_print_results();
//
//	////plot_wait_on_join(state);
//	////plot_shutdown(state);
//
//	return 0;
//}

snow_main();
