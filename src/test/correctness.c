/*
 * This file includes tests to ensure the correctness of the library and
 * the computations it performs.
 *
 * Common acronyms you will encounter:
 * 	ITEM -- Imaginary Time Evolution Method
 * 	HO 	-- Harmonic Oscilllator
 * 	SCIM -- Self-Consistency Iteration Method
 * 	FDM -- Finite Difference Method
 *
 * Currently the following is being tested (in the order they appear in
 * this source file):
 * 	[x] bucket array
 * 	[x] priority queue
 * 	[ ] complex hermitian bandmat vector multiplication
 * 	[ ] HO function sampling
 * 	[x] quadgk 1D numerical integration
 * 	[x] quadgk_vec 1D numerical integration
 * 	[ ] ------ 2D numerical integration
 * 	[ ] ------ 3D numerical integration
 * 	[ ] FDM vs HO hamiltonian solving
 * 		[ ] particle in a box
 * 		[ ] perturbed HO potential
 * 	[ ] ITEM vs SCIM groundstate finding
 */

#include <sbmf/sbmf.h>
#include <sbmf/methods/quadgk_vec.h>
#include <sbmf/methods/quadgk_vec_inl.h>
#include <sbmf/methods/find_groundstate.h>
#include <sbmf/math/functions.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/math/matrix.h>
#include <sbmf/math/find_eigenpairs.h>
#include <sbmf/math/grid.h>
#include <sbmf/memory/prioqueue.h>
#include <sbmf/memory/bucketarray.h>
#include <sbmf/debug/profile.h>

#include <omp.h>
#include <plot/plot.h>
#define SNOW_ENABLED
#include "snow.h"

#include <math.h>
#include <stdio.h>

#define DISABLE_SLOW_TESTS 0

describe(random) {
	it ("works") {
		const u32 N = 5;
		for (u32 i = 0; i < N; ++i) {
			for (u32 j = i; j < N; ++j) {
				log_info("%u,%u", i,j);
			}
		}
		log_info("next one");
		for (u32 i = 0; i < (N*(N+1))/2; ++i) {
			log_info("%u", i);
		}

		//
		//					0,0		0,1		0,2
		//					--		1,1		1,2
		//					--		--		2,2
		//
		//
		//					u32 i = (H.size-1)*(H.size-(c-r)) + r;
		//
		//
		//					(0) 0,0		(1)0,1		(3)0,2
		//					--			(2)1,1		(4)1,2
		//					--			--			(5)2,2
		//
		//					total number of unique combinations in a n*n grid
		//					is given by the arithmetic series
		//						1 + 2 + 3 + ... + n = n(n+1)/2
		//
		//					i  |  r,c
		//					---+-----
		//					0  |  0,0
		//					1  |  0,1
		//					2  |  0,2
		//					3  |  1,1
		//					4  |  1,2
		//					5  |  2,2
		//
		//
		//					i = (n-1)*(n - (c-r)) + r
		//
		//					0,0 -> 2*(3-(0-0))+0 = 6
		//					0,1 -> 2*(3-(1-0))+0 = 4
		//					0,2 -> 2*(3-(2-0))+0 = 2
		//					1,1 -> 2*(3-(1-1))+1 = 7
		//					1,2 -> 2*(3-(2-1))+1 = 5
		//					2,2 -> 2*(3-(2-2))+2 = 8
		//

	}
}

/* bucket array */
describe (bucketarray) {
	before_each() { sbmf_init(); }
	after_each() { sbmf_shutdown(); }
	it ("creating, settings and getting values")  {
		struct barray* ba = barray_new(2, sizeof(i32));
		asserteq(ba->bucket_count, 1);

		i32 i = 0;

		i = 10; barray_set(ba, 0, &i);
		i =  9; barray_set(ba, 1, &i);
		asserteq(ba->bucket_count, 1);

		i =  8; barray_set(ba, 2, &i);
		i =  7; barray_set(ba, 3, &i);
		asserteq(ba->bucket_count, 2);

		i =  6; barray_set(ba, 6, &i);
		asserteq(ba->bucket_count, 4);

		asserteq(*(i32*)barray_get(ba,0), 10);
		asserteq(*(i32*)barray_get(ba,1),  9);
		asserteq(*(i32*)barray_get(ba,2),  8);
		asserteq(*(i32*)barray_get(ba,3),  7);
		asserteq(*(i32*)barray_get(ba,6),  6);
	}
}

/* priority queue */

static bool pqcmp_des(void* a, void* b) {
	return *((i32*)a) < *((i32*)b);
}

static bool pqcmp_asc(void* a, void* b) {
	return *((i32*)a) > *((i32*)b);
}

describe (pq) {
	before_each() { sbmf_init(); }
	after_each()  { sbmf_shutdown(); }

	it ("pushes and pops (des)") {
		struct prioqueue* pq = prioqueue_new(4, sizeof(i32), pqcmp_des);
		for (i32 i = 0; i < 100; ++i) {
			prioqueue_push(pq, &i);
		}

		i32 val;
		for (i32 i = 0; i < 100; ++i) {
			prioqueue_pop(pq, &val);
			asserteq(val, i);
		}

		asserteq(pq->mem->bucket_count, 100/4);
	}
	it ("pushes and pops (asc)") {
		struct prioqueue* pq = prioqueue_new(4, sizeof(i32), pqcmp_asc);
		for (i32 i = 0; i < 1000; ++i) {
			prioqueue_push(pq, &i);
			prioqueue_push(pq, &i);
			prioqueue_pop(pq, NULL);
		}

		i32 val;
		for (i32 i = 999; i >= 0; --i) {
			prioqueue_pop(pq, &val);
			asserteq(val, i);
		}

		//asserteq(pq->mem->bucket_count, 100/4);
	}
}

/* complex hermitian bandmat vector multiplication */
describe (complex_hermitian_bandmat) {
	before_each() { sbmf_init(); }
	after_each()  { sbmf_shutdown(); }

	it ("identity vector mult.") {
		const u32 N = 10;

		struct complex_hermitian_bandmat identity = complex_hermitian_bandmat_new(1,N);
		for (u32 r = 0; r < identity.size; ++r) {
			u32 i = complex_hermitian_bandmat_index(identity, r,r);
			identity.data[i] = 1.0;
		}

		c64 in[N];
		c64 out[N];
		for (u32 i = 0; i < N; ++i)
			in[i] = i;

		complex_hermitian_bandmat_mulv(out, identity, in);

		for (u32 i = 0; i < N; ++i) {
			asserteq(out[i], in[i]);
		}
	}

	it ("more complicated mat. vector mult.") {
		const u32 N = 5;

		struct complex_hermitian_bandmat bm = complex_hermitian_bandmat_new_zero(2,N);
		COMPLEX_HERMITIAN_BANDMAT_FOREACH(bm, r,c) {
			u32 i = complex_hermitian_bandmat_index(bm, r,c);
			bm.data[i] = (f64)r + (f64)((i32)r-(i32)c)*I;
		}

		c64 in[N];
		c64 out[N];
		c64 expected[N];
		expected[0] = 0-1*I;
		expected[1] = 3-2*I;
		expected[2] = 11-2*I;
		expected[3] = 25-2*I;
		expected[4] = 25+3*I;

		for (u32 i = 0; i < N; ++i)
			in[i] = (f64)i;

		complex_hermitian_bandmat_mulv(out, bm, in);

		for (u32 i = 0; i < N; ++i) {
			asserteq(out[i], expected[i]);
		}
	}

	it ("full mat. vector mult.") {
		const u32 N = 3;

		struct complex_hermitian_bandmat bm = complex_hermitian_bandmat_new_zero(N,N);
		COMPLEX_HERMITIAN_BANDMAT_FOREACH(bm, r,c) {
			u32 i = complex_hermitian_bandmat_index(bm, r,c);
			bm.data[i] = (f64)r + (f64)((i32)r-(i32)c)*I;
		}

		/*
		 *			 0		-I		-2I
		 *			+I		1		1-I
		 *			+2I		1+I		2
		 *
		 *
		 *			-I-4I		-5I
		 *			1+2-2I = 	3-2I
		 *			1+I+4		5+I
		 *
		 */

		c64 in[N];
		c64 out[N];
		c64 expected[N];
		expected[0] = 0-5*I;
		expected[1] = 3-2*I;
		expected[2] = 5+I;

		for (u32 i = 0; i < N; ++i)
			in[i] = (f64)i;

		complex_hermitian_bandmat_mulv(out, bm, in);

		for (u32 i = 0; i < N; ++i) {
			asserteq(out[i], expected[i]);
		}
	}
}

/* HO function sampling */
describe (hofunctionsampling) {
}

/* quadgk 1D numerical integration */

static void check_quadgk_converge(integration_result res, f64 expected) {
	bool correct_ans = f64_compare(res.integral, expected, 1e-9);
	if (!res.converged || !correct_ans) {
		printf("Integral failed to converge or got wrong answer:\n\tconverged: %d\n\tintegral: %lf\n\terror: %lf\n\texpected: %lf\n\tevals: %d\n", res.converged, res.integral, res.error, expected, res.performed_evals);
	}

asserteq(correct_ans && res.converged, true);
}

/* quadgk_vec 1D numerical integration */

void x2_vec(f64* out, f64* in, u32 len, void* p) {
	SBMF_UNUSED(p);
	for (u32 i = 0; i < len; ++i)
		out[i] = in[i]*in[i];
}
void expx_vec(f64* out, f64* in, u32 len, void* p) {
	SBMF_UNUSED(p);
	for (u32 i = 0; i < len; ++i)
		out[i] = exp(in[i]);
}
void expnx_vec(f64* out, f64* in, u32 len, void* p) {
	SBMF_UNUSED(p);
	for (u32 i = 0; i < len; ++i)
		out[i] = exp(-in[i]);
}
void expnabsx_vec(f64* out, f64* in, u32 len, void* p) {
	SBMF_UNUSED(p);
	for (u32 i = 0; i < len; ++i)
		out[i] = exp(-fabs(in[i]));
}
void sinx_vec(f64* out, f64* in, u32 len, void* p) {
	SBMF_UNUSED(p);
	for (u32 i = 0; i < len; ++i)
		out[i] = sin(in[i]);
}

describe (quad_gk_vec_numerical_integration){
	before_each() { sbmf_init(); }
	after_each() { sbmf_shutdown(); }

	integration_settings settings = {
		.gk = gk7,
		.abs_error_tol = 1e-10,
		.rel_error_tol = 1e-10,
		.max_evals = 500,
	};

	it ("x2, 0 -> 2") {
		integration_result res = quadgk_vec(x2_vec, 0, 2, settings);
		check_quadgk_converge(res, 8.0/3.0);
	}

	it ("x2, 2 -> 0") {
		integration_result res = quadgk_vec(x2_vec, 2, 0, settings);
		check_quadgk_converge(res, -8.0/3.0);
	}

	it ("x2, -2 -> 0") {
		integration_result res = quadgk_vec(x2_vec, -2, 0, settings);
		check_quadgk_converge(res, 8.0/3.0);
	}

	it ("x2, 0 -> -2") {
		integration_result res = quadgk_vec(x2_vec, 0, -2, settings);
		check_quadgk_converge(res, -8.0/3.0);
	}

	it ("x2, -2 -> 2") {
		integration_result res = quadgk_vec(x2_vec, -2, 2, settings);
		check_quadgk_converge(res, 2*8.0/3.0);
	}

	it ("expnx, 0 -> inf") {
		integration_result res = quadgk_vec(expnx_vec, 0, INFINITY, settings);
		check_quadgk_converge(res, 1.0);
	}

	it ("expnx, inf -> 0") {
		integration_result res = quadgk_vec(expnx_vec, INFINITY, 0, settings);
		check_quadgk_converge(res, -1.0);
	}

	it ("expx, -inf -> 0") {
		integration_result res = quadgk_vec(expx_vec, -INFINITY, 0, settings);
		check_quadgk_converge(res, 1.0);
	}

	it ("expx, 0 -> -inf") {
		integration_result res = quadgk_vec(expx_vec, 0, -INFINITY, settings);
		check_quadgk_converge(res, -1.0);
	}

	it ("|expnx|, -inf -> inf") {
		integration_result res = quadgk_vec(expnabsx_vec, -INFINITY, INFINITY, settings);
		check_quadgk_converge(res, 2.0);
	}

	it ("|expnx|, inf -> -inf") {
		integration_result res = quadgk_vec(expnabsx_vec, INFINITY, -INFINITY, settings);
		check_quadgk_converge(res, -2.0);
	}

	it ("sin, 0 -> 2pi") {
		integration_result res = quadgk_vec(sinx_vec, 0, 2*M_PI, settings);
		check_quadgk_converge(res, 0.0);
	}

	it ("sin, 0 -> pi") {
		integration_result res = quadgk_vec(sinx_vec, 0, M_PI, settings);
		check_quadgk_converge(res, 2.0);
	}
}

/* quadgk_vec_inl 1D numerical integration */

describe (quad_gk_vec_inl_numerical_integration){
	before_each() { sbmf_init(); }
	after_each() { sbmf_shutdown(); }

	integration_settings settings = {
		.gk = gk7,
		.abs_error_tol = 1e-10,
		.rel_error_tol = 1e-10,
		.max_evals = 500,
	};

	it ("x2, 0 -> 2") {
		integration_result res = quadgk_vec_inl(x2_vec, 0, 2, settings);
		check_quadgk_converge(res, 8.0/3.0);
	}

	it ("x2, 2 -> 0") {
		integration_result res = quadgk_vec_inl(x2_vec, 2, 0, settings);
		check_quadgk_converge(res, -8.0/3.0);
	}

	it ("x2, -2 -> 0") {
		integration_result res = quadgk_vec_inl(x2_vec, -2, 0, settings);
		check_quadgk_converge(res, 8.0/3.0);
	}

	it ("x2, 0 -> -2") {
		integration_result res = quadgk_vec_inl(x2_vec, 0, -2, settings);
		check_quadgk_converge(res, -8.0/3.0);
	}

	it ("x2, -2 -> 2") {
		integration_result res = quadgk_vec_inl(x2_vec, -2, 2, settings);
		check_quadgk_converge(res, 2*8.0/3.0);
	}

	it ("expnx, 0 -> inf") {
		integration_result res = quadgk_vec_inl(expnx_vec, 0, INFINITY, settings);
		check_quadgk_converge(res, 1.0);
	}

	it ("expnx, inf -> 0") {
		integration_result res = quadgk_vec_inl(expnx_vec, INFINITY, 0, settings);
		check_quadgk_converge(res, -1.0);
	}

	it ("expx, -inf -> 0") {
		integration_result res = quadgk_vec_inl(expx_vec, -INFINITY, 0, settings);
		check_quadgk_converge(res, 1.0);
	}

	it ("expx, 0 -> -inf") {
		integration_result res = quadgk_vec_inl(expx_vec, 0, -INFINITY, settings);
		check_quadgk_converge(res, -1.0);
	}

	it ("|expnx|, -inf -> inf") {
		integration_result res = quadgk_vec_inl(expnabsx_vec, -INFINITY, INFINITY, settings);
		check_quadgk_converge(res, 2.0);
	}

	it ("|expnx|, inf -> -inf") {
		integration_result res = quadgk_vec_inl(expnabsx_vec, INFINITY, -INFINITY, settings);
		check_quadgk_converge(res, -2.0);
	}

	it ("sin, 0 -> 2pi") {
		integration_result res = quadgk_vec_inl(sinx_vec, 0, 2*M_PI, settings);
		check_quadgk_converge(res, 0.0);
	}

	it ("sin, 0 -> pi") {
		integration_result res = quadgk_vec_inl(sinx_vec, 0, M_PI, settings);
		check_quadgk_converge(res, 2.0);
	}
}

/* ------ 2D numerical integration */
/* ------ 3D numerical integration */

/* FDM vs HOB particle in a box */
describe (fdm_vs_hob_particle_in_a_box) {
	before_each() { sbmf_init(); }
	after_each() { sbmf_shutdown(); }
}

/* FDM vs HOB perturbed harmonic oscillator */

static f64 ho_perturbed_potential(f64* x, i32 n, void* data) {
	SBMF_UNUSED(n);
	SBMF_UNUSED(data);
	return ho_potential(x,1,0) + gaussian(*x,0,0.2);
}

typedef f64 pot_func(f64*,i32,void*);
struct integrand_params {
	u32 n[2];
	pot_func* pot;
};

typedef void pot_func_vec(f64*,f64*,u32,void*);
struct integrand_params_vec {
	u32 n[2];
	pot_func_vec* pot;
};

/* ITEM vs SCIM groundstate finding */

c64 guess(f64* v, i32 n) {
	SBMF_UNUSED(n);
	return gaussian(*v, 0, 0.2);
}

void guess_vec(c64* out, f64* in, u32 len) {
	for (u32 i = 0; i < len; ++i) {
		out[i] = gaussian(in[i],0,0.2);
	}
}

f64 linear_hamiltonian_pot(f64* v, i32 n, c64 u) {
	SBMF_UNUSED(u);
	return ho_perturbed_potential(v, n, NULL);
}

void linear_hamiltonian_vec_pot(f64* out, f64* in_x, c64* in_u, u32 len) {
	SBMF_UNUSED(in_u);
	for (u32 i = 0; i < len; ++i) {
		out[i] = gaussian(in_x[i],0,0.2);
	}
}

f64 non_linear_hamiltonian_pot(f64* v, i32 n, c64 u) {
	/* assuming 1d */
	return ho_potential(v, n, 0) + cabs(u)*cabs(u);
}

void non_linear_hamiltonian_vec_pot(f64* out, f64* in_x, c64* in_u, u32 len) {
	SBMF_UNUSED(in_x);
	for (u32 i = 0; i < len; ++i) {
		out[i] = cabs(in_u[i])*cabs(in_u[i]);
	}
}

void debug_callback(struct scim_settings settings, c64* wavefunction) {
	plot_init(800, 600, "scim debug");
	const u32 N = 128;
	const f64 L = 5.0;
	f32 pdata[N];
	f32 wdata[N];
	f32 rewdata[N];
	f32 imwdata[N];
	sample_space sp = make_linspace(1, -L/2.0, L/2.0, N);

	for (u32 i = 0; i < N; ++i) {
		c64 sample = hob_sample(wavefunction, settings.num_basis_functions, sp.points[i]);

		f64 pdataout;
		f64 x = sp.points[i];
		non_linear_hamiltonian_vec_pot(&pdataout, &x, &sample, 1);
		pdata[i] = (f32)pdataout;

		c64 c = cabs(sample);
		wdata[i] = c*c;
		rewdata[i] = creal(sample);
		imwdata[i] = cimag(sample);
	}
	push_line_plot(&(plot_push_desc){
			.space = &sp,
			.data = pdata,
			.label = "potential",
			});
	push_line_plot(&(plot_push_desc){
			.space = &sp,
			.data = wdata,
			.label = "abs",
			});

	push_line_plot(&(plot_push_desc){
			.space = &sp,
			.data = rewdata,
			.label = "re",
			});
	push_line_plot(&(plot_push_desc){
			.space = &sp,
			.data = imwdata,
			.label = "im",
			});
	plot_update_until_closed();
	plot_shutdown();
}

describe(item_vs_scim_groundstate_finding) {
	before_each() { sbmf_init(); }
	after_each() { sbmf_shutdown(); }

	it ("perturbed HO potential") {
		const f64 L = 10.0;
		const u32 N = 256;
		struct grid space = generate_grid(1,
				(f64[]){-L/2.0},
				(f64[]){+L/2.0},
				(i32[]){N});
		struct item_settings item_settings = {
			.g = space,
			.max_iterations = 1e7,
			.error_tol = 1e-10,
			.dt = 0.001,
		};
		struct scim_settings scim_settings = {
			.num_basis_functions = 32,
			.max_iterations = 1e7,
			.error_tol = 1e-10,
		};

		struct gss_result item_res = item(item_settings, linear_hamiltonian_pot, guess);
		log_info("\nitem:\niterations: %d\nerror: %e", item_res.iterations, item_res.error);

		struct gss_result hob_res = ho_scim(scim_settings, linear_hamiltonian_vec_pot, guess_vec);
		log_info("\nhob:\niterations: %d\nerror: %e", hob_res.iterations, hob_res.error);
#if 0
		{
			plot_init(800, 600, "fdm groundstate");
			f32 pdata[N];
			sample_space sp = make_linspace(1, -L/2.0, L/2.0, N);

			for (u32 i = 0; i < N; ++i) {
				f64 x = sp.points[i];
				pdata[i] = (f32) ho_perturbed_potential(&x, 1, NULL);
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = pdata,
					.label = "potential",
					});

			for (u32 i = 0; i < N; ++i) {
				c64 c = cabs(item_res.wavefunction[i]);
				pdata[i] = c*c;
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = pdata,
					.label = "item groundstate",
					});

			for (u32 i = 0; i < N; ++i) {
				c64 c = cabs(hob_sample(hob_res.wavefunction, scim_settings.num_basis_functions, sp.points[i]));
				pdata[i] = c*c;
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = pdata,
					.label = "hob groundstate",
					});

			plot_update_until_closed();
			plot_shutdown();
		}
#endif
	}

	it ("non-linear hamiltonian") {
		const f64 L = 10.0;
		const u32 N = 512;
		struct grid space = generate_grid(1,
				(f64[]){-L/2.0},
				(f64[]){+L/2.0},
				(i32[]){N});
		struct item_settings item_settings = {
			.g = space,
			.max_iterations = 1e9,
			.error_tol = 1e-9,
			.dt = 0.0001,
		};
		struct scim_settings scim_settings = {
			.num_basis_functions = 32,
			.max_iterations = 1e9,
			.error_tol = 1e-7,
			//.measure_every = 40,
			//.dbgcallback = debug_callback,
		};

		PROFILE_BEGIN("entire item");
		struct gss_result item_res = item(item_settings, non_linear_hamiltonian_pot, guess);
		PROFILE_END("entire item");
		log_info("\nitem:\niterations: %d\nerror: %e", item_res.iterations, item_res.error);

		PROFILE_BEGIN("entire ho_scim");
		struct gss_result hob_res = ho_scim(scim_settings, non_linear_hamiltonian_vec_pot, guess_vec);
		PROFILE_END("entire ho_scim");
		log_info("\nhob:\niterations: %d\nerror: %e", hob_res.iterations, hob_res.error);

#if 0
		{
			plot_init(800, 600, "fdm groundstate");
			f32 pdata[N];
			sample_space sp = make_linspace(1, -L/2.0, L/2.0, N);

			/*
			for (u32 i = 0; i < N; ++i) {
				f64 x = sp.points[i];
				pdata[i] = (f32) ho_perturbed_potential(&x, 1, NULL);
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = pdata,
					.label = "potential",
					});
			*/

			for (u32 i = 0; i < N; ++i) {
				c64 c = cabs(item_res.wavefunction[i]);
				pdata[i] = c*c;
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = pdata,
					.label = "item groundstate",
					});

			c64 sample_out[N];
			f64 sample_in[N];
			for (u32 i = 0; i < N; ++i) {
				sample_in[i] = (f64) sp.points[i];
			}
			hob_sample_vec(hob_res.wavefunction, scim_settings.num_basis_functions, sample_out, sample_in, N);
			for (u32 i = 0; i < N; ++i) {
				f64 c = cabs(sample_out[i]);
				pdata[i] = c*c;
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = pdata,
					.label = "hob groundstate",
					});

			plot_update_until_closed();
			plot_shutdown();
		}
#endif
	}
}

snow_main();
