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
 * 	[ ] 2 component GP solving
 */

#include <sbmf/sbmf.h>
#include <sbmf/methods/quadgk_vec.h>
#include <sbmf/methods/item.h>
#include <sbmf/methods/nlse_solver.h>
#include <sbmf/methods/grosspitaevskii.h>
//#include <sbmf/methods/gp2c_gsl.h>
#include <sbmf/methods/best_meanfield.h>
#include <sbmf/math/functions.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/math/matrix.h>
#include <sbmf/math/find_eigenpairs.h>
#include <sbmf/math/grid.h>
#include <sbmf/memory/prioqueue.h>
#include <sbmf/memory/bucketarray.h>
#include <sbmf/debug/profile.h>
#include <sbmf/math/manybodystate.h>

#include <omp.h>
#include <plot/plot.h>
#define SNOW_ENABLED
#include "snow.h"

#include <math.h>
#include <stdio.h>
#include <errno.h>

#define DISABLE_SLOW_TESTS 0

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


/* quadgk 1D numerical integration */
static void check_quadgk_converge(integration_result res, f64 expected) {
	bool correct_ans = f64_compare(res.integral, expected, 1e-9);
	if (!res.converged || !correct_ans) {
		printf("Integral failed to converge or got wrong answer:\n\tconverged: %d\n\tintegral: %lf\n\terror: %e\n\texpected: %lf\n\tevals: %d\n", res.converged, res.integral, res.error, expected, res.performed_evals);
	}

	asserteq(correct_ans && res.converged, true);
}

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
		.abs_error_tol = 1e-7,
		.rel_error_tol = 1e-7,
		.max_evals = 1e5,
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

/* ITEM vs SCIM groundstate finding */
c64 item_guess(f64* v, i32 n) {
	SBMF_UNUSED(n);
	return gaussian(*v, 0, 0.2);
}

f64 linear_hamiltonian_pot(f64* v, i32 n, c64 u) {
	SBMF_UNUSED(u);
	return ho_potential(v,1,0) + gaussian(*v,0,0.2);
}

void linear_hamiltonian_vec_pot(f64* out, f64* in_x, f64* in_u, u32 len) {
	SBMF_UNUSED(in_u);
	for (u32 i = 0; i < len; ++i) {
		out[i] = gaussian(in_x[i],0,0.2);
	}
}

f64 non_linear_hamiltonian_pot(f64* v, i32 n, c64 u) {
	/* assuming 1d */
	return ho_potential(v,1,0) + gaussian(*v, 0, 0.2)  - 3.0*cabs(u)*cabs(u);
}

void non_linear_hamiltonian_vec_pot(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count], void* userdata) {
	for (u32 i = 0; i < len; ++i) {
		out[i] = gaussian(in_x[i],0,0.2) - 3.0*fabs(in_u[i])*fabs(in_u[i]);
	}
}

describe(item_vs_gp) {
	before_each() { sbmf_init(); }
	after_each() { sbmf_shutdown(); }

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
			.dt = 1e-4,
		};

		struct nlse_settings nlse_settings = {
			.num_basis_funcs = 32,
			.max_iterations = 1e9,
			.error_tol = 1e-9,
			.gk = gk15,
			.basis = ho_basis,
			.zero_threshold=1e-10,
		};

		struct gss_result item_res = item(item_settings, non_linear_hamiltonian_pot, item_guess);

		struct nlse_result nlse_res = nlse_solver(nlse_settings, 1, &(struct nlse_component) {
					.op = non_linear_hamiltonian_vec_pot,
				});

		f64 nlse_sample[N];
		nlse_settings.basis.sample(nlse_res.coeff_count, nlse_res.coeff, N, nlse_sample, space.points);

		for (u32 i = 0; i < N; ++i) {
			f64 c1 = fabs(nlse_sample[i]);
			f64 c2 = cabs(item_res.wavefunction[i]);
			asserteq(fabs(c1*c1 - c2*c2) < 1e-2, true);
		}

#if 0
		{
			plot_init(800, 600, "fdm groundstate");
			f32 pdata[N];
			sample_space sp = make_linspace(1, -L/2.0, L/2.0, N);

			for (u32 i = 0; i < N; ++i) {
				c64 c = cabs(item_res.wavefunction[i]);
				pdata[i] = c*c;
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = pdata,
					.label = "item groundstate",
					});

			f64 sample_out[N];
			f64 sample_in[N];
			for (u32 i = 0; i < N; ++i) {
				sample_in[i] = (f64) sp.points[i];
			}
			gp2c_settings.basis.sample(gp2c_res.coeff_count, gp2c_res.coeff, N, sample_out, sample_in);
			for (u32 i = 0; i < N; ++i) {
				f64 c = fabs(sample_out[i]);
				pdata[i] = c*c;
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = pdata,
					.label = "gp2c",
					});

			plot_update_until_closed();
			plot_shutdown();
		}
#endif
	}
}






















#define GAA_2COMP (-4.0)
#define GBB_2COMP (-4.0)
#define GAB_2COMP (+4.0)

void op_2comp_a(const u32 len, f64 out[static len],
		f64 in_x[static len], const u32 component_count,
		f64 in_u[static len*component_count], void* userdata) {
	for (u32 i = 0; i < len; ++i) {
		f64 ca = fabs(in_u[i]);
		f64 cb = fabs(in_u[len + i]);
		out[i] =
			+ GAA_2COMP*ca*ca
            + GAB_2COMP*cb*cb;
	}
}

void op_2comp_b(const u32 len, f64 out[static len],
		f64 in_x[static len], const u32 component_count,
		f64 in_u[static len*component_count], void* userdata) {
	for (u32 i = 0; i < len; ++i) {
		f64 ca = fabs(in_u[i]);
		f64 cb = fabs(in_u[len + i]);
		out[i] =
			+ GBB_2COMP*cb*cb
            + GAB_2COMP*ca*ca;
	}
}

void perturbation(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count], void* userdata) {
    assert(component_count == 0);
    for (u32 i = 0; i < len; ++i) {
        out[i] = gaussian(in_x[i], 0, 0.2);
    }
}


//describe(2comp) {
//	before_each() { sbmf_init(); }
//	after_each()  { sbmf_shutdown(); }
//
//	it ("?") {
//		struct gp2c_settings gp2c_settings = {
//			.num_basis_functions = 16,
//			.max_iterations = 1e9,
//			.error_tol = 1e-9,
//			.gk = gk15,
//			.basis = ho_basis,
//			.ho_potential_perturbation = perturbation,
//		};
//
//		struct gp2c_component comps[2] = {
//			[0] = {
//				.op = op_2comp_a,
//			},
//			[1] = {
//				.op = op_2comp_b,
//			},
//		};
//
//		struct gp2c_result gp2c_res = gp2c(gp2c_settings, 2, comps);
//
//#if 1
//		{
//			const f64 L = 10.0;
//			const u32 N = 512;
//
//			plot_init(800, 600, "2comp");
//			f32 pdata[N];
//			sample_space sp = make_linspace(1, -L/2.0, L/2.0, N);
//
//			for (u32 i = 0; i < N; ++i) {
//				f64 x = sp.points[i];
//				pdata[i] = (f32) ho_potential(&x,1,0) + gaussian(x,0,0.2);
//			}
//			push_line_plot(&(plot_push_desc){
//					.space = &sp,
//					.data = pdata,
//					.label = "potential",
//					});
//
//			f64 sample_in[N];
//			for (u32 i = 0; i < N; ++i) {
//				sample_in[i] = (f64) sp.points[i];
//			}
//
//			for (u32 i = 0; i < gp2c_res.component_count; ++i) {
//				f64 sample_out[N];
//				gp2c_settings.basis.sample(gp2c_res.coeff_count, &gp2c_res.coeff[i*gp2c_res.coeff_count], N, sample_out, sample_in);
//
//				for (u32 i = 0; i < N; ++i) {
//					f64 c = fabs(sample_out[i]);
//					pdata[i] = c*c;
//				}
//				push_line_plot(&(plot_push_desc){
//						.space = &sp,
//						.data = pdata,
//						.label = plot_snprintf("comp %u", i),
//						.offset = gp2c_res.energy[i],
//						});
//			}
//
//
//			plot_update_until_closed();
//			plot_shutdown();
//		}
//#endif
//	}
//}
//
//
//
//
//





































#if 0
describe (2comp_scim) {
	before_each(){sbmf_init();}
	after_each(){sbmf_shutdown();}

	it ("reproduces 1 component solutions") {
		struct gp2c_settings settings = {
			.num_basis_functions = 8,
			.max_iterations = 1e7,
			.error_tol = 1e-8,
		};
		struct gp2c_result res = gp2c(settings, gp2c_op_a, gp2c_op_b);

#if 0
		const u32 states_to_include = 5;
		/* First order correction of the many-body wavefunction */
		{
			struct eigen_result eres_a = find_eigenpairs_sparse(res.hamiltonian_a, states_to_include, EV_SMALLEST_RE);

			struct Vijkl_params params = {
				.coeff_count = settings.num_basis_functions,
			};

			integration_settings int_settings = {
				.gk = gk7,
				.abs_error_tol = 1e-10,
				.rel_error_tol = 1e-10,
				.max_evals = 1e7,
				.userdata = &params,
			};

			/*
			 * |0> = |N,0,...>
			 * |mn> = |N-2,0,...,0,1,0,...,0,1,0,...>
			 * 					   ^-_ m:th  ^-_ n:th
			 * |> = c0|...> + c1|...> + ...
			 * <x|> = c0<x|...> + c1<x|...> + ...
			 * <x|abc...> = phi_a(x)phi_b(x)phi_c(x)...
			 */
			struct manybody_state {
				struct {
					u32 particles_in_state;
					u32 state_index;
				} indices[3];
				u32 index_count;
			};

			/*
			 *
			 *		x x x x
			 *		x x x x
			 *		x x x x
			 *		x x x x
			 *		1 + 2 + 3 + 4 = 10
			 */
			const u32 state_count = 1 + (states_to_include-1)*(states_to_include-1+1)/2;
			u32 index = 0;
			struct manybody_state states[state_count];
			f64 coeffs[state_count];

			states[0] = (struct manybody_state) {
				.indices[0] = {
					.particles_in_state = PARTICLE_COUNT,
					.state_index = 0
				},
				.index_count = 1,
			};
			coeffs[0] = 1;
			index++;

			/* We know the inner product of single subst. and the groundstate is zero.
			 * Only have to account for double substitutions
			 */

			/* m,n here denotes the location of the excited orbitals */
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 n = m; n < states_to_include; ++n) {
					/* Currently disregarding double counting. */

					params.i = m;

					params.j = n;
					params.coeff_i = &eres_a.eigenvectors[m * eres_a.num_eigenpairs];
					params.coeff_j = &eres_a.eigenvectors[n * eres_a.num_eigenpairs];
					params.coeff_0 = &eres_a.eigenvectors[0];

					integration_result int_res = quadgk_vec(Vijkl_integrand, -INFINITY, INFINITY, int_settings);

					assert(int_res.converged);
					f64 Vmn00 = int_res.integral;

					f64 inner_product = GAA * sqrt(PARTICLE_COUNT*(PARTICLE_COUNT-1)) * Vmn00;

					printf("%u -- (%u,%u)\n", index,m,n);
					if (m == n) {
						states[index] = (struct manybody_state) {
							.indices[0] = {
								.particles_in_state = PARTICLE_COUNT-2,
								.state_index = 0
							},
							.indices[1] = {
								.particles_in_state = 2,
								.state_index = m,
							},
							.index_count = 2,
						};

						f64 energy_diff =
							PARTICLE_COUNT*eres_a.eigenvalues[0]
							- (PARTICLE_COUNT-2)*eres_a.eigenvalues[0]
							- eres_a.eigenvalues[m]
							- eres_a.eigenvalues[n];
						coeffs[index] = sqrt(2)*inner_product/energy_diff;
					} else {
						states[index] = (struct manybody_state) {
							.indices[0] = {
								.particles_in_state = PARTICLE_COUNT-2,
								.state_index = 0
							},
							.indices[1] = {
								.particles_in_state = 1,
								.state_index = m,
							},
							.indices[2] = {
								.particles_in_state = 1,
								.state_index = n,
							},
							.index_count = 3,
						};
						f64 energy_diff =
							PARTICLE_COUNT*eres_a.eigenvalues[0]
							- (PARTICLE_COUNT-2)*eres_a.eigenvalues[0]
							- eres_a.eigenvalues[m]
							- eres_a.eigenvalues[n];
						coeffs[index] = inner_product/energy_diff;
					}

					index++;
				}
			}

			f64 sum = 0.0;
			for (u32 i = 0; i < state_count; ++i) {
				sum += coeffs[i]*coeffs[i];
			}
			f64 scaling = 1.0/sum;
			for (u32 i = 0; i < state_count; ++i) {
				coeffs[i] *= scaling;
			}

			for (u32 i = 0; i < state_count; ++i) {
				printf("%.2lf -- ", coeffs[i]);
				for (u32 j = 0; j < states[i].index_count; ++j) {
					printf("(%u,%u)",
							states[i].indices[j].state_index,
							states[i].indices[j].particles_in_state
							);
				}
				printf("\n");
			}

			printf("-------------------------------\n");

			const u32 N = 256;
			plot_init(800, 600, "gp2c");
			f32 potdata[N], adata[N];
			memset(adata, 0, sizeof(adata[0])*N);
			sample_space sp = make_linspace(1, -5, 5.0, N);

			for (u32 i = 0; i < N; ++i) {
				f64 x = sp.points[i];
				potdata[i] = (f32) ho_perturbed_potential(&x, 1, NULL);
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = potdata,
					.label = "potential",
					});

			f64 sample_in[N];
			for (u32 i = 0; i < N; ++i) {
				sample_in[i] = (f64) sp.points[i];
			}

			c64 out[N];
			for (u32 i = 0; i < N; ++i) {
				out[i] = 0;
			}
			for (u32 i = 0; i < state_count; ++i) {

				c64 mp_out[N];
				for (u32 j = 0; j < N; ++j) {
					mp_out[j] = 1.0;
				}

				for (u32 j = 0; j < states[i].index_count; ++j) {
					c64 sample_out[N];
					hob_sample_vec(
							&eres_a.eigenvectors[states[i].indices[j].state_index * eres_a.num_eigenpairs],
							settings.num_basis_functions,
							sample_out,
							sample_in,
							N);

					for (u32 k = 0; k < N; ++k) {
						mp_out[k] *= cpow(sample_out[k], states[i].indices[j].particles_in_state);
					}
				}

				for (u32 j = 0; j < N; ++j) {
					out[j] += coeffs[i]*mp_out[j];
				}
			}

			for (u32 i = 0; i < N; ++i) {
				f32 tmp = (f32)fabs(out[i]);
				adata[i] = tmp*tmp;
			}

			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = adata,
					.label = "a",
					.offset = res.energy_a,
					});

			plot_update_until_closed();
			plot_shutdown();




		}
#endif

#if 1
		{
			const u32 N = 256;
			plot_init(800, 600, "gp2c");
			f32 potdata[N], adata[N], bdata[N];
			sample_space sp = make_linspace(1, -5, 5.0, N);

			for (u32 i = 0; i < N; ++i) {
				f64 x = sp.points[i];
				potdata[i] = (f32) ho_perturbed_potential(&x, 1, NULL);
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = potdata,
					.label = "potential",
					});

			c64 sample_out_a[N];
			c64 sample_out_b[N];
			f64 sample_in[N];
			for (u32 i = 0; i < N; ++i) {
				sample_in[i] = (f64) sp.points[i];
			}
			hob_sample_vec(res.coeff_a, settings.num_basis_functions, sample_out_a, sample_in, N);
			hob_sample_vec(res.coeff_b, settings.num_basis_functions, sample_out_b, sample_in, N);
			for (u32 i = 0; i < N; ++i) {
				f64 ca = fabs(sample_out_a[i]);
				f64 cb = fabs(sample_out_b[i]);
				adata[i] = ca*ca;
				bdata[i] = cb*cb;
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = adata,
					.label = "a",
					.offset = res.energy_a,
					});
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = bdata,
					.label = "b",
					.offset = res.energy_b,
					});

			plot_update_until_closed();
			plot_shutdown();
		}
#endif
	}
}
#endif

























void bestmf_perturbation(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count]) {
	assert(component_count == 0);
	for (u32 i = 0; i < len; ++i) {
		out[i] = gaussian(in_x[i],0,0.2);
	}
}

static u32 bestmf_particle_count = 100;
static f64 bestmf_interaction_strength = 0.0;

void bestmf_gp2c_op_a(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count]) {
	assert(component_count == 1);
	f64 gaa = bestmf_interaction_strength/(bestmf_particle_count - 1);

#pragma omp simd
	for (u32 i = 0; i < len; ++i) {
		f64 ca = fabs(in_u[i]);
		out[i] = gaa*(bestmf_particle_count-1)*ca*ca;
	}
}

//describe (bestmf) {
//	before_each(){sbmf_init();}
//	after_each(){sbmf_shutdown();}
//
//	it ("?") {
//		struct gp2c_settings settings = {
//			.num_basis_functions = 8,
//			.max_iterations = 1e7,
//			.error_tol = 1e-8,
//			.dbgcallback = bestmf_debug_callback,
//			.measure_every = 0,
//			.ho_potential_perturbation = bestmf_perturbation,
//		};
//		struct gp2c_component component = {
//			.op = bestmf_gp2c_op_a,
//		};
//		struct gp2c_result res = gp2c(settings, 1, &component);
//
//		struct eigen_result eres_a = find_eigenpairs_sparse(res.hamiltonian[0], 2, EV_SMALLEST_RE);
//		c64_normalize(&eres_a.eigenvectors[0], &eres_a.eigenvectors[0], settings.num_basis_functions);
//		c64_normalize(&eres_a.eigenvectors[settings.num_basis_functions], &eres_a.eigenvectors[settings.num_basis_functions], settings.num_basis_functions);
//#if 0
//		{
//			const u32 N = 256;
//			plot_init(800, 600, "bestmf gp2c");
//			f32 potdata[N], adata[N], b0data[N], b1data[N];
//			sample_space sp = make_linspace(1, -5, 5.0, N);
//
//			for (u32 i = 0; i < N; ++i) {
//				f64 x = sp.points[i];
//				potdata[i] = (f32) ho_perturbed_potential(&x, 1, NULL);
//			}
//			push_line_plot(&(plot_push_desc){
//					.space = &sp,
//					.data = potdata,
//					.label = "potential",
//					});
//
//			c64 sample_out_a[N];
//			c64 sample_out_b0[N];
//			c64 sample_out_b1[N];
//			f64 sample_in[N];
//			for (u32 i = 0; i < N; ++i) {
//				sample_in[i] = (f64) sp.points[i];
//			}
//			hob_sample_vec(res.coeff, settings.num_basis_functions, sample_out_a, sample_in, N);
//			hob_sample_vec(&eres_a.eigenvectors[0], settings.num_basis_functions, sample_out_b0, sample_in, N);
//			hob_sample_vec(&eres_a.eigenvectors[settings.num_basis_functions], settings.num_basis_functions, sample_out_b1, sample_in, N);
//			for (u32 i = 0; i < N; ++i) {
//				f64 ca = cabs(sample_out_a[i]);
//				adata[i] = ca*ca;
//
//				f64 cb0 = cabs(sample_out_b0[i]);
//				b0data[i] = cb0*cb0;
//
//				f64 cb1 = cabs(sample_out_b1[i]);
//				b1data[i] = cb1*cb1;
//			}
//			push_line_plot(&(plot_push_desc){
//					.space = &sp,
//					.data = adata,
//					.label = "a",
//					.offset = res.energy[0],
//					});
//
//			push_line_plot(&(plot_push_desc){
//					.space = &sp,
//					.data = b0data,
//					.label = "b0",
//					.offset = eres_a.eigenvalues[0]
//					});
//
//			push_line_plot(&(plot_push_desc){
//					.space = &sp,
//					.data = b1data,
//					.label = "b1",
//					.offset = eres_a.eigenvalues[1]
//					});
//
//			plot_update_until_closed();
//			plot_shutdown();
//		}
//#endif
//
//		sbmf_log_info("---------------------------------");
//
//
//		sbmf_log_info("0 -- energy: %lf", eres_a.eigenvalues[0]);
//		for (u32 i = 0; i < settings.num_basis_functions; ++i) {
//			printf("%.1lf\t",cabs(eres_a.eigenvectors[i]));
//		}
//		printf("\n");
//		sbmf_log_info("1 -- energy: %lf", eres_a.eigenvalues[1]);
//		for (u32 i = 0; i < settings.num_basis_functions; ++i) {
//			printf("%.1lf\t",cabs(eres_a.eigenvectors[eres_a.num_eigenpairs + i]));
//		}
//		printf("\n");
//
//		u32 particle_count = 100;
//		f64 g = (-2.5)/(particle_count-1);
//		find_best_meanfield_occupations(particle_count, g, settings.num_basis_functions,
//				&eres_a.eigenvectors[0],
//				&eres_a.eigenvectors[settings.num_basis_functions],
//				eres_a.eigenvalues[0],
//				eres_a.eigenvalues[1]
//				);
//	}
//}

//describe(bestmf_fig_1) {
//	before_each(){sbmf_init();}
//	after_each(){sbmf_shutdown();}
//
//	it ("?") {
//		struct gp2c_settings settings = {
//			.num_basis_functions = 32,
//			.max_iterations = 1e8,
//			.error_tol = 1e-8,
//			.dbgcallback = bestmf_debug_callback,
//			.measure_every = 0,
//			.ho_potential_perturbation = bestmf_perturbation,
//		};
//		struct gp2c_component component = {
//			.op = bestmf_gp2c_op_a,
//		};
//		bestmf_particle_count = 500;
//
//		const f64 gvals[] = {
//			(-4.0)/(bestmf_particle_count-1),
//			(-3.5)/(bestmf_particle_count-1),
//			(-2.5)/(bestmf_particle_count-1),
//			(-1.5)/(bestmf_particle_count-1),
//			(-1.0)/(bestmf_particle_count-1),
//			(-0.5)/(bestmf_particle_count-1),
//			( 0.5)/(bestmf_particle_count-1)
//		};
//
//		FILE* bestfile = fopen("bestmf_data_solved", "w");
//		for (u32 i = 0; i < ARRLEN(gvals); ++i) {
//			static char buf[50];
//			snprintf(buf, 50, "output/bestmf_data_%.1lf", gvals[i]*(bestmf_particle_count-1));
//			FILE* datafile = fopen(buf, "w");
//			if (!datafile) {
//				sbmf_log_error("Unable to open log file: %s", buf);
//				sbmf_log_error("errno: (%d) %s", errno, strerror(errno));
//			}
//
//			bestmf_interaction_strength = gvals[i];
//
//			struct gp2c_result res = gp2c(settings, 1, &component);
//			struct eigen_result eres_a = find_eigenpairs_sparse(res.hamiltonian[0], 2, EV_SMALLEST_RE);
//			c64_normalize(&eres_a.eigenvectors[0], &eres_a.eigenvectors[0], settings.num_basis_functions);
//			c64_normalize(&eres_a.eigenvectors[settings.num_basis_functions], &eres_a.eigenvectors[settings.num_basis_functions], settings.num_basis_functions);
//
//			for (u32 j = 0; j <= bestmf_particle_count; j += 5) {
//				const f64 energy = best_meanfield_energy(settings.num_basis_functions,
//											&eres_a.eigenvectors[0],
//											&eres_a.eigenvectors[settings.num_basis_functions],
//											j,
//											bestmf_particle_count - j,
//											gvals[i]);
//				fprintf(datafile, "%lf\t%lf\n", (f64)j/(f64)bestmf_particle_count, energy/(f64)bestmf_particle_count);
//			}
//			fclose(datafile);
//
//			struct best_meanfield_results bmfres = find_best_meanfield_occupations(bestmf_particle_count, gvals[i], settings.num_basis_functions,
//												&eres_a.eigenvectors[0],
//												&eres_a.eigenvectors[settings.num_basis_functions],
//												eres_a.eigenvalues[0],
//												eres_a.eigenvalues[1]);
//			fprintf(bestfile, "%lf\t%lf\n", (f64)bmfres.occupation_1/(f64)bestmf_particle_count, bmfres.energy/(f64)bestmf_particle_count);
//		}
//		fclose(bestfile);
//	}
//}


describe(check_full_energy_calc) {
	before_each(){sbmf_init();}
	after_each(){sbmf_shutdown();}

	it ("does it") {
		struct gp_settings settings = {
			.num_basis_funcs = 16,
			.basis = ho_basis,
			.component_count = 1,
			.occupations = (u32[]){2},
			.g0 = (f64[]){0.0},
		};

		struct nlse_result res = grosspitaevskii(settings);
		asserteq(f64_compare(res.energy[0], 0.5, 1e-5), true);
		asserteq(f64_compare(full_energy(settings, res), 1.0, 1e-5), true);
	}
}
















snow_main();
