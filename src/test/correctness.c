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
 * 	[ ] priority queue
 * 	[x] quadgk 1D numerical integration
 * 	[x] quadgk_vec 1D numerical integration
 * 	[ ] ------ 2D numerical integration
 * 	[ ] ------ 3D numerical integration
 * 	[ ] FDM vs HO hamiltonian solving
 * 		[ ] particle in a box
 * 		[ ] perturbed HO potential
 * 	[ ] ITEM vs SCIM groundstate finding
 */

#define SNOW_ENABLED
#include "snow.h"

#include <math.h>

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

#include <plot/plot.h>

#define DISABLE_SLOW_TESTS 0

/* priority queue */

static bool pqcmp(void* a, void* b) {
	return *((i32*)a) < *((i32*)b);
}

describe (pq) {
	before_each() { sbmf_init(); }
	after_each() { sbmf_shutdown(); }

	it ("pushes and pops") {
		prioqueue* pq = prioqueue_new(10, sizeof(i32), pqcmp);
		i32 i = 0;

		i = 10; prioqueue_push(pq, &i);
		i =  9; prioqueue_push(pq, &i);
		i =  8; prioqueue_push(pq, &i);
		i =  7; prioqueue_push(pq, &i);
		i =  6; prioqueue_push(pq, &i);
		i =  5; prioqueue_push(pq, &i);
		i =  4; prioqueue_push(pq, &i);
		i =  3; prioqueue_push(pq, &i);
		i =  2; prioqueue_push(pq, &i);
		i =  1; prioqueue_push(pq, &i);

		asserteq(*(i32*)prioqueue_top(pq), 1 ); prioqueue_pop(pq);
		asserteq(*(i32*)prioqueue_top(pq), 2 ); prioqueue_pop(pq);
		asserteq(*(i32*)prioqueue_top(pq), 3 ); prioqueue_pop(pq);
		asserteq(*(i32*)prioqueue_top(pq), 4 ); prioqueue_pop(pq);
		asserteq(*(i32*)prioqueue_top(pq), 5 ); prioqueue_pop(pq);
		asserteq(*(i32*)prioqueue_top(pq), 6 ); prioqueue_pop(pq);
		asserteq(*(i32*)prioqueue_top(pq), 7 ); prioqueue_pop(pq);
		asserteq(*(i32*)prioqueue_top(pq), 8 ); prioqueue_pop(pq);
		asserteq(*(i32*)prioqueue_top(pq), 9 ); prioqueue_pop(pq);
		asserteq(*(i32*)prioqueue_top(pq), 10); prioqueue_pop(pq);
	}
}

/* quadgk 1D numerical integration */

static void check_quadgk_converge(integration_result res, f64 expected) {
	bool correct_ans = f64_compare(res.integral, expected, 1e-9);
	if (!res.converged || !correct_ans) {
		printf("Integral failed to converge or got wrong answer:\n\tconverged: %d\n\tintegral: %lf\n\terror: %lf\n\texpected: %lf\n", res.converged, res.integral, res.error, expected);
	}

	asserteq(correct_ans && res.converged, true);
}

/* quadgk_vec 1D numerical integration */

void x2_vec(f64* out, f64* in, u32 len, void* p) {
	for (u32 i = 0; i < len; ++i)
		out[i] = in[i]*in[i];
}
void expx_vec(f64* out, f64* in, u32 len, void* p) {
	for (u32 i = 0; i < len; ++i)
		out[i] = exp(in[i]);
}
void expnx_vec(f64* out, f64* in, u32 len, void* p) {
	for (u32 i = 0; i < len; ++i)
		out[i] = exp(-in[i]);
}
void expnabsx_vec(f64* out, f64* in, u32 len, void* p) {
	for (u32 i = 0; i < len; ++i)
		out[i] = exp(-fabs(in[i]));
}
void sinx_vec(f64* out, f64* in, u32 len, void* p) {
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
	return ho_potential(x,1,0) + gaussian(*x,0,0.2);
}

static void ho_perturbed_potential_vec(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i) {
		out[i] = ho_potential(&in[i],1,0) + gaussian(in[i],0,0.2);
	}
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

static f64 hob_integrand(f64 x, void* data) {
	struct integrand_params* params = data;
	return
		ho_eigenfunction((i32[]){params->n[0]}, &x, 1) *
		params->pot(&x,1,NULL) *
		ho_eigenfunction((i32[]){params->n[1]}, &x, 1);
}

static void hob_integrand_vec(f64* out, f64* in, u32 len, void* data) {
	struct integrand_params_vec* params = data;

	f64 eig1[len];
	f64 eig2[len];
	f64 pot[len];

	ho_eigenfunction_vec(params->n[0], eig1, in, len);
	ho_eigenfunction_vec(params->n[1], eig2, in, len);
	params->pot(pot, in, len, data);

	for (u32 i = 0; i < len; ++i) {
		out[i] = eig1[i]*eig2[i]*pot[i];
	}
}

static struct eigen_result fdm_solve(pot_func* pot, f64 L, f64 N) {
	struct grid space = generate_grid(1,
			(f64[]){-L/2.0},
			(f64[]){+L/2.0},
			(i32[]){N});
	hermitian_bandmat mat = construct_finite_diff_mat(space.pointcounts[0], space.dimensions, space.deltas);

	// Invert and scale matrix to kinetic energy term
	for (mat_size_t i = 0; i < mat.base.rows*mat.base.cols; ++i) {
		mat.base.data[i] = -0.5 * mat.base.data[i];
	}

	// Loop through main diagonal and add potential energy term
	for (mat_size_t i = 0; i < mat.size; ++i) {
		mat_size_t index = mat.size*(mat.bandcount-1) + i;
		mat.base.data[index] += pot(&space.points[i], space.dimensions, NULL);
	}

	asserteq(mat_is_valid(mat.base), true);

	// Solve eigenvalue problem for hamiltonian
	return find_eigenpairs_sparse(mat, 3, EV_SMALLEST_RE);
}

static struct eigen_result hob_solve(pot_func_vec* pot, f64 N) {
	hermitian_bandmat T = construct_ho_kinetic_matrix(N);

	struct integrand_params_vec params;
	integration_settings settings = {
		.gk = gk7,
		.abs_error_tol = 1e-10,
		.rel_error_tol = 1e-10,
		.max_evals = 1e4,
		.userdata = &params
	};

	params.pot = pot;
	for (u32 r = 0; r < T.size; ++r) {
		for (u32 c = r; c < T.size; ++c) {
			params.n[0] = r;
			params.n[1] = c;
			integration_result res = quadgk_vec(hob_integrand_vec, -INFINITY, INFINITY, settings);
			asserteq(res.converged, true);

			u32 i = (T.size-1)*(T.size-(c-r)) + r;
			T.base.data[i] += res.integral;
		}
	}

	asserteq(mat_is_valid(T.base), true);

	// Solve eigenvalue problem for hamiltonian and return
	return find_eigenpairs_sparse(T, 3, EV_SMALLEST_RE);
}

describe (fdm_vs_hob) {
	before_each() { sbmf_init(); }
	after_each() { sbmf_shutdown(); }

	it ("particle in a box") {
	}

	it ("perturbed harmonic oscillator") {
		const f64 L = 5.0;
		const i32 N = 64;
		struct eigen_result fdm_res = fdm_solve(ho_perturbed_potential, L, N);
		struct eigen_result hob_res = hob_solve(ho_perturbed_potential_vec, N/2);
	}
}

/* ITEM vs SCIM groundstate finding */

c64 guess(f64* v, i32 n) {
	return gaussian(*v, 0, 0.2);
}


void guess_vec(c64* out, f64* in, u32 len) {
	for (u32 i = 0; i < len; ++i) {
		out[i] = gaussian(in[i],0,0.2);
	}
}

f64 linear_hamiltonian_pot(f64* v, i32 n, c64 u) {
	return ho_perturbed_potential(v, n, NULL);
}

void linear_hamiltonian_vec_pot(f64* out, f64* in_x, c64* in_u, u32 len) {
	for (u32 i = 0; i < len; ++i) {
		out[i] = ho_perturbed_potential(&in_x[i], 1, NULL);
	}
}

f64 non_linear_hamiltonian_pot(f64* v, i32 n, c64 u) {
	/* assuming 1d */
	return ho_potential(v, n, 0) + 3*cabs(u)*cabs(u);
}

void non_linear_hamiltonian_vec_pot(f64* out, f64* in_x, c64* in_u, u32 len) {
	ho_potential_vec(out, in_x, len);
	for (u32 i = 0; i < len; ++i) {
		out[i] += 3*cabs(in_u[i])*cabs(in_u[i]);
	}
}

describe(item_vs_scim_groundstate_finding) {
	before_each() { sbmf_init(); }
	after_each() { sbmf_shutdown(); }

	it ("perturbed HO potential") {
		const f64 L = 5.0;
		const i32 N = 128;
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
		c64_normalize(item_res.wavefunction, space.total_pointcount);

		struct gss_result hob_res = scim(scim_settings, linear_hamiltonian_vec_pot, guess_vec);
		log_info("\nhob:\niterations: %d\nerror: %e", hob_res.iterations, hob_res.error);
#if 1
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
				pdata[i] = 30*c*c;
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
		const i32 N = 256;
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
			.num_basis_functions = 64,
			.max_iterations = 1e7,
			.error_tol = 1e-10,
		};

		struct gss_result item_res = item(item_settings, non_linear_hamiltonian_pot, guess);
		log_info("\nitem:\niterations: %d\nerror: %e", item_res.iterations, item_res.error);
		c64_normalize(item_res.wavefunction, space.total_pointcount);

		struct gss_result hob_res = scim(scim_settings, non_linear_hamiltonian_vec_pot, guess_vec);
		log_info("\nhob:\niterations: %d\nerror: %e", hob_res.iterations, hob_res.error);

#if 1
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
				pdata[i] = 10*c*c;
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
