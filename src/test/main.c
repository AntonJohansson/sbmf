#include <sbmf/groundstate_solver/groundstate_solver.h>
#include <sbmf/common/profile.h>
#include <sbmf/common/fdm.h>
#include <sbmf/common/grid.h>
#include <sbmf/common/eigenproblem.h>

#include <stdio.h>
#include <assert.h>
#include <unistd.h>

#include "plotting/plt.h"

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

static u64 factorial(u64 n) {
	assert(n <= 20); // largest factorial supported by u64

	if (n < 1) {
		return 1;
	} else {
		return n*factorial(n-1);
	}
}

static inline f64 hermite_poly(i32 n, f64 x) {
	// H_n(z) = (-1)^n * exp(z^2) * (d^n/dz^n) (exp(-z^2))
	//
	// Explicit form:
	// 		H_n(x) = n! sum_(m=0)^(floor(n/2)) (-1)^m/(m!*(n-2m)!) * (2x)^(n-2m)

	f64 sum = 0.0;
	i32 upper_bound = n/2;
	for (i32 m = 0; m <= upper_bound; ++m) {
		sum += pow(-1,m) / (factorial(m) * factorial(n-2*m)) * pow(2*x, n-2*m);
	}

	return factorial(n)*sum;
}

static inline f64 psi_ho(i32 states[], f64 point[], i32 dims) {
	f64 sum = 0.0;
	for (i32 i = 0; i < dims; ++i)
		sum += exp(-point[i]*point[i]/2.0) * hermite_poly(states[i], point[i]);
	return sum;
}

static inline f64 energy_ho(i32 states[], i32 dims) {
	f64 sum = 0.0;
	for (i32 i = 0; i < dims; ++i) 
		sum += states[i];
	
	return sum + dims/2.0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline void test_eigenvalue_solving_harmonic_osc_1d() {
	printf("-- Running 1D harmonic osc test\n");

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
			fdm.bands[idx] += (c64) harmonic_osc_potential(&g.points[i], g.dimensions, 0.0);
		}	
	}
	PROFILE_END("gen. fdm  ho 1d");


	PROFILE_BEGIN("e.v. prob. ho 1d");
	f64* eigenvalues = malloc(fdm.size * sizeof(f64));
	c64* eigenvectors = malloc(fdm.order * sizeof(c64));
	eigvp_bandmat(eigenvalues, eigenvectors, fdm);
	PROFILE_END("e.v. prob. ho 1d");

	PlotState* state = plt_init();
	f32 x[g.total_pointcount];
	f32 v[g.total_pointcount];
	for (i32 i = 0; i < g.total_pointcount; ++i) {
		x[i] = g.points[i];
		v[i] = harmonic_osc_potential(&g.points[i], g.dimensions, 0.0);
	}
	plt_1d(state, x, v, g.total_pointcount);

	f64 expected_answer[g.total_pointcount];
	for (i32 i = 0; i < 10; ++i) {
		normalize_function_c64(&eigenvectors[i*g.total_pointcount], g.total_pointcount);
		
		for (i32 j = 0; j < g.total_pointcount; ++j) {
			//expected_answer[j] = psi_ho_1d(i, g.points[j]);
			expected_answer[j] = psi_ho(&i, &g.points[j], g.dimensions);
		}
		normalize_function_f64(expected_answer, g.total_pointcount);

		for (i32 j = 0; j < g.total_pointcount; ++j) {
			f64 expected = expected_answer[j];
			f64 got = creal(eigenvectors[i*g.total_pointcount + j]);
			if (fabs(expected-got) > 0.01) {
				fprintf(stderr, "\t- Eigenvector calculation failed: got %lf; expected: %lf\n", got, expected);
			}
		}

		f64 expected_energy = energy_ho(&i, g.dimensions);
		if (fabs(eigenvalues[i] - expected_energy) > 0.01) {
			fprintf(stderr, "\t - Eigenvalue calculation failed: got %lf, expected %lf\n", eigenvalues[i], expected_energy);
		}

		f32 y[g.total_pointcount];
		f32 z[g.total_pointcount];
		for (i32 j = 0; j < g.total_pointcount; ++j) {
			//y[j] = E_ho_1d(i) + expected_answer[j];
			y[j] = energy_ho(&i, g.dimensions) + expected_answer[j];
			z[j] = eigenvalues[i] + creal(eigenvectors[i*g.total_pointcount + j]);
		}

		plt_1d(state, x, y, g.total_pointcount);
		plt_1d(state, x, z, g.total_pointcount);
	}

	plt_wait_on_join(state);
	plt_shutdown(state);

	(void)normalize_function_f32;
	free_grid(g);
}

static inline void test_eigenvalue_solving_harmonic_osc_2d() {
	printf("-- Running 2D harmonic osc test\n");

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
		v[i] = harmonic_osc_potential(&g.points[2*i], g.dimensions, 0.0);
		u[i] = -10*creal(eigenvectors[1*g.total_pointcount + i]);
		e[i] = psi_ho((i32[]){1,0}, &g.points[2*i], g.dimensions);
	}

	for (i32 i = 0; i < g.total_pointcount; ++i) {
		if (fabs(e[i] - u[i]) > 0.01) {
			fprintf(stderr, "\t- Eigenvector calculation %d failed: got %lf; expected %lf\n", i, e[i], u[i]);
		}
	}

	normalize_function_f32(&e[0], g.total_pointcount);
	for (i32 i = 0; i < g.total_pointcount; ++i) {
		e[i] = 10*e[i];
	}

	PlotState* state = plt_init();
	plt_2d(state, x, y, v, g.total_pointcount);
	plt_2d(state, x, y, u, g.total_pointcount);
	plt_2d(state, x, y, e, g.total_pointcount);
	plt_wait_on_join(state);
	plt_shutdown(state);
	
	
	
	
	

	
	
	
	

	free(eigenvalues);
	free(eigenvectors);
	free_grid(g);
}

static inline void test_eigenvalue_solving_harmonic_osc_3d() { }
static inline void test_eigenvalue_solving_particle_in_a_box_1d() { }
static inline void test_eigenvalue_solving_particle_in_a_box_2d() { }
static inline void test_eigenvalue_solving_particle_in_a_box_3d() { }

static inline void test_eigenvalue_solving_periodic_pot_1d() {
	printf("-- Running 1D periodic test\n");

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


	PlotState* state = plt_init();

	f32 x[g.total_pointcount];
	f32 v[g.total_pointcount];
	f32 u[g.total_pointcount];
	for (i32 i = 0; i < g.total_pointcount; ++i) {
		x[i] = g.points[i];
		v[i] = periodic_pot(&g.points[i], g.dimensions, 0.0);
	}

	plt_1d(state, x, v, g.total_pointcount);

	for (i32 i = 0; i < 5; ++i) {
		normalize_function_c64(&eigenvectors[i*g.total_pointcount], g.total_pointcount);
		for (i32 j = 0; j < g.total_pointcount; ++j) {
			u[j] = 10*creal(eigenvectors[i*g.total_pointcount + j]);
		}

		plt_1d(state, x, u, g.total_pointcount);
	}


	plt_wait_on_join(state);
	plt_shutdown(state);

	(void)normalize_function_f32;
	free_grid(g);
}

static inline void test_eigenvalue_solving() {
	printf("-- General band matrix\n");

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

	// produces by matlab, warning COLUMN MAJOR
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

	// check answer
	for (int i = 0; i < bm.size; ++i) {
		if (fabs(eigvals[i] - eigvals_answer[i]) > 0.01) {
			fprintf(stderr, "\t- Eigenvalue calculation %d failed: got %lf; expected %lf\n", i, eigvals[i], eigvals_answer[i]);
		}
	}

	(void)eigvecs_answer;
	// Check eigenvectors
	for (int i = 0; i < bm.size; ++i) {
		for (int j = 0; j < bm.size; ++j) {
			c64 got = eigvecs[i*bm.size + j]; 
			c64 expected = eigvecs_answer[j*bm.size + i];
			if (cabs(got)-cabs(expected) > 0.01) {
				fprintf(stderr, "\t- Eigenvector calculation failed: got %lf + %lf i; expected: %lf + %lf i\n", creal(got), cimag(got), creal(expected), cimag(expected));
			}
		}
	}
}

int main(int argc, char** argv) {
	assert(hermite_poly(0,0) == 1);
	assert(hermite_poly(1,5) == 2*5);
	assert(hermite_poly(2,5) == 4*5*5-2);

	assert(factorial(0) == 1);
	assert(factorial(5) == 120);
	assert(factorial(10) == 3628800);
	assert(factorial(20) == 2432902008176640000);

	printf("Testing eigenvalue solving\n"),
	//test_eigenvalue_solving();
	test_eigenvalue_solving_harmonic_osc_1d();
	test_eigenvalue_solving_periodic_pot_1d();
	//test_eigenvalue_solving_harmonic_osc_2d();

	//PlotState* state = plt_init();

	//const i32 N = 128;

	//grid g = generate_grid(1, 
	//		(f64[]){-5, -5},
	//		(f64[]){ 5,  5},
	//		 (i32[]){ N,  	 	 N     }
	//		);

	//gss_settings settings = {
	//	.g = g,
	//	.max_iterations = 2000,
	//	.error_tol = 1e-10,
	//	.dt = 0.01,

	//	.measure_every = 0,
	//};

	//PROFILE_BEGIN("item");
	//gss_result res = item_execute(settings,
	//															potential,
	//															initial_guess);
	//PROFILE_END("item");
	//printf(
	//		"ITEM results:\n"
	//		"\titerations: %d/%d\n"
	//		"\terror: %e\n"
	//		"\terror tol: %e\n",
	//	res.iterations,
	//	res.settings.max_iterations,
	//	res.error,
	//	res.settings.error_tol
	//);


	//float x[g.total_pointcount];
	//float y[g.total_pointcount];
	//float v[g.total_pointcount];
	//for (i32 i = 0; i < g.total_pointcount; ++i) {
	//	x[i] = (float)g.points[i];
	//	y[i] = (float)creal(res.wavefunction[i]);
	//	v[i] = (float)potential(&g.points[i], 1, res.wavefunction[i]);
	//}

	//plt_1d(state, x, v, g.total_pointcount); 
	//plt_1d(state, x, y, g.total_pointcount); 


	//PROFILE_BEGIN("gen. fdm");
	//bandmat fdm;
	//{
	//	// NOTE: g.pointcounts[0] forces it to be square!
	//	fdm = generate_fd_matrix(g.pointcounts[0], pow(2,g.dimensions), g.dimensions, g.deltas);

	//	for (i32 i = 0; i < fdm.size*fdm.bandcount; ++i) {
	//		fdm.bands[i] = -0.5*fdm.bands[i];
	//	}

	//	for (i32 i = 0; i < fdm.size; ++i) {
	//		i32 idx = fdm.size*(fdm.bandcount-1) + i;
	//		fdm.bands[idx] += (c64) potential(&g.points[i], 1, res.wavefunction[i]);
	//	}	
	//}
	//PROFILE_END("gen. fdm");


	//PROFILE_BEGIN("e.v. prob.");
	//{
	//	f64* eigenvalues = malloc(fdm.size * sizeof(f64));
	//	c64* eigenvectors = malloc(fdm.order * sizeof(c64));

	//	eigvp_bandmat(eigenvalues, eigenvectors, fdm);

	//	float x[fdm.size];
	//	float u[fdm.size];
	//	for (int i = 1; i < 6; ++i) {
	//		for (int j = 0; j < fdm.size; ++j) {
	//			x[j] = (float)g.points[j];
	//			u[j] = (float)eigenvalues[i] +  3*(float)creal(eigenvectors[j*fdm.size + i]);
	//		}
	//		plt_1d(state, x,u, fdm.size);
	//	}
	//}
	//PROFILE_END("e.v. prob.");




	//gss_free_result(res);
	//free_grid(g);

	profile_print_results();

	////plt_wait_on_join(state);
	////plt_shutdown(state);

	return 0;
}
