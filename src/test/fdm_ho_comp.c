#include <math.h>
#include <sbmf/common/matrix.h>
#include <sbmf/common/eigenproblem.h>
#include <sbmf/basis/harmonic_oscillator.h>
#include <sbmf/numerical_integration/quadgk.h>
#include <plot/plot.h>

#define PLOT_HO_HO_PERT 0

void c64_normalize(c64* data, u32 size) {
	f64 sum = 0.0;
	for (u32 i = 0; i < size; ++i) {
		f64 absval = cabs(data[i]);
		sum += absval*absval;
	}
	for (u32 i = 0; i < size; ++i) {
		data[i] *= (1.0/sqrt(sum));
	}
}

static inline f64 gaussian(f64 x, f64 mu, f64 sigma) {
	return 1.0/(sigma*sqrt(M_2_PI)) * exp(-x*x/(2*sigma*sigma));
}

static inline f64 ho_perturbed_potential(f64* x, i32 n, void* data) {
	return ho_potential(x,1,0) + gaussian(*x,0,0.2);
}

typedef f64 pot_func(f64*,i32,void*);
typedef struct {
	u32 n[2];
	pot_func* pot;
} hob_integrand_params;

static inline f64 hob_integrand(f64 x, void* data) {
	hob_integrand_params* params = data;
	return
		ho_eigenfunction((i32[]){params->n[0]}, &x, 1) *
		params->pot(&x,1,NULL) *
		ho_eigenfunction((i32[]){params->n[1]}, &x, 1);
}

static eig_result fdm_solve(pot_func* pot, f64 L, f64 N) {
	grid space = generate_grid(1,
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
	eig_result res = eig_sparse_bandmat(mat, 3, EV_SMALLEST_RE);
	return res;
}

static eig_result hob_solve(pot_func* pot, f64 N) {
	hermitian_bandmat T = construct_ho_kinetic_matrix(N);

	hob_integrand_params params;
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
			integration_result res = quadgk(hob_integrand, -INFINITY, INFINITY, settings);
			asserteq(res.converged, true);

			u32 i = (T.size-1)*(T.size-(c-r)) + r;
			T.base.data[i] += res.integral;
		}
	}

	asserteq(mat_is_valid(T.base), true);

	// Solve eigenvalue problem for hamiltonian
	eig_result res = eig_sparse_bandmat(T, 3, EV_SMALLEST_RE);
	return res;
}

describe(fdm_ho_comp) {
	before_each() {
		sbmf_init();
	}
	after_each() {
		sbmf_shutdown();
	}
	it ("fdm") {
	}

	it ("ho fdm comp") {
		const f64 L = 5.0;
		const i32 N = 128;
		eig_result fdm_res = fdm_solve(ho_perturbed_potential, L, N);
		eig_result hob_res = hob_solve(ho_perturbed_potential, N/2);

		{
			plot_init(800, 600, "ho fdm comp");
			f32 pdata[N];
			c64 cdata[N];
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

			for (u32 i = 0; i < fdm_res.num_eigenpairs; ++i) {
				c64_normalize(&fdm_res.eigenvectors[i*fdm_res.points_per_eigenvector], fdm_res.points_per_eigenvector);
				for (u32 j = 0; j < sp.pointcount; ++j) {
					c64 c = fdm_res.eigenvectors[i*fdm_res.points_per_eigenvector + j];
					pdata[j] = cabs(c)*cabs(c);
				}
				push_line_plot(&(plot_push_desc){
						.space = &sp,
						.data = pdata,
						.label = plot_snprintf("fdm %d", i),
						.offset = cabs(fdm_res.eigenvalues[i]),
						});
			}
			for (u32 i = 0; i < hob_res.num_eigenpairs; ++i) {
				memset(cdata, 0, sizeof(cdata));
				for (u32 j = 0; j < hob_res.points_per_eigenvector; ++j) {
					c64 coeff = hob_res.eigenvectors[i*hob_res.points_per_eigenvector + j];

					// Get the j:th basis function
					for (u32 k = 0; k < sp.pointcount; ++k) {
						f64 x = sp.points[k];
						cdata[k] += coeff*ho_eigenfunction((i32[]){j}, &x, 1);
					}
				}
				c64_normalize(cdata, sp.pointcount);

				for (u32 j = 0; j < sp.pointcount; ++j) {
					pdata[j] = cabs(cdata[j])*cabs(cdata[j]);
				}

				push_line_plot(&(plot_push_desc){
						.space = &sp,
						.data = pdata,
						.label = plot_snprintf("hob %d", i),
						.offset = cabs(hob_res.eigenvalues[i]),
						});
			}
			plot_update_until_closed();
			plot_shutdown();
		}
	}
}
