#include <math.h>
#include <sbmf/common/matrix.h>
#include <sbmf/common/eigenproblem.h>
#include <sbmf/basis/harmonic_oscillator.h>
#include <plot/plot.h>

void f64_normalize(f64* data, u32 size) {
	f64 sum = 0.0;
	for (u32 i = 0; i < size; ++i) {
		f64 absval = fabs(data[i]);
		sum += absval*absval;
	}
	for (u32 i = 0; i < size; ++i) {
		data[i] *= (1.0/sqrt(sum));
	}
}
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

f64 exact_pob_eigval(f64 L, u32 n) {
	return (n*n*M_PI*M_PI)/(2*L*L);
}

f64 exact_pob_eigfunc(f64 x, f64 L, u32 n) {
	if (n % 2 == 0)
		return sqrt(2.0/L)*sin((n*M_PI/L)*x);
	else
		return sqrt(2.0/L)*cos((n*M_PI/L)*x);
}

f64 potential_well(f64 x, f64 L) {
	if (fabs(x) < L/2.0)
		return 0.0;
	else
		return 10.0;
}

f64 potential_well_smooth(f64 x, f64 L) {
	return 100*pow(x/(L/2.0), 32);
}

static inline f64 pob_integrand(f64 x, void* data) {
	u32* params = data;
	return
		ho_eigenfunction((i32[]){params[0]}, &x, 1) *
		potential_well(x, 1.0) *
		ho_eigenfunction((i32[]){params[1]}, &x, 1);
}

static inline f64 pob_integrand_smooth(f64 x, void* data) {
	u32* params = data;
	return
		ho_eigenfunction((i32[]){params[0]}, &x, 1) *
		potential_well_smooth(x, 1.0) *
		ho_eigenfunction((i32[]){params[1]}, &x, 1);
}

static inline f64 ho_integrand(f64 x, void* data) {
	u32* params = data;
	return
		ho_eigenfunction((i32[]){params[0]}, &x, 1) *
		ho_potential(&x,1,0)*
		ho_eigenfunction((i32[]){params[1]}, &x, 1);
}

describe(particle_in_a_box_1D) {
	const f64 L = 1.0;
	const u32 max_n = 5;
	const u32 sample_points = 128;
	f64 exact_pob_func[max_n][sample_points];
	f64 exact_pob_energy[max_n];

	grid space;

	before_each() {
		sbmf_init();
		space = generate_grid(1,
				(f64[]){-L/2.0},
				(f64[]){+L/2.0},
				(i32[]){sample_points});
		for (u32 n = 0; n < max_n; ++n) {
			exact_pob_energy[n] = exact_pob_eigval(L, n);
			for (u32 i = 0; i < space.total_pointcount; ++i) {
				f64 x = space.points[i];
				exact_pob_func[n][i] = exact_pob_eigfunc(x, L, n+1);
			}
			f64_normalize(&exact_pob_func[n][0], space.total_pointcount);
		}
	}

	after_each() {
		sbmf_shutdown();
	}

	it ("FDM hamiltonian") {
		hermitian_bandmat mat = construct_finite_diff_mat(space.pointcounts[0], space.dimensions, space.deltas);

		// Invert and scale matrix to kinetic energy term
		for (mat_size_t i = 0; i < mat.base.rows*mat.base.cols; ++i) {
			mat.base.data[i] = -0.5 * mat.base.data[i];
		}

		// Loop through main diagonal and add potential energy term
		for (mat_size_t i = 0; i < mat.size; ++i) {
			mat_size_t index = mat.size*(mat.bandcount-1) + i;
			mat.base.data[index] += potential_well(space.points[i], L);
		}

		asserteq(mat_is_valid(mat.base), true);

		// Solve eigenvalue problem for hamiltonian
		eig_result res = eig_sparse_bandmat(mat, max_n, EV_SMALLEST_RE);

		// Check results
		for (u32 n = 0; n < res.num_eigenpairs; ++n) {
			c64_normalize(&res.eigenvectors[n*res.points_per_eigenvector], res.points_per_eigenvector);
			for (u32 i = 0; i < res.points_per_eigenvector; ++i) {
				f64 diff = cabs(res.eigenvectors[n*res.points_per_eigenvector + i]) - fabs(exact_pob_func[n][i]);
				asserteq(fabs(diff) < 0.05, true);
			}
		}
	}

	it ("HO basis hamiltonian -- ho") {
		hermitian_bandmat T = construct_ho_kinetic_matrix(max_n+2);

		u32 params[2];
		integration_settings settings = {
			.gk = gk7,
			.abs_error_tol = 1e-10,
			.rel_error_tol = 1e-10,
			.max_evals = 1e4,
			.userdata = params
		};

		for (u32 r = 0; r < T.size; ++r) {
			for (u32 c = r; c < T.size; ++c) {
				params[0] = r;
				params[1] = c;
				integration_result res = quadgk(ho_integrand, -INFINITY, INFINITY, settings);
				asserteq(res.converged, true);

				u32 i = (T.size-1)*(T.size-(c-r)) + r;
				T.base.data[i] += res.integral;
			}
		}

		asserteq(mat_is_valid(T.base), true);

		// Solve eigenvalue problem for hamiltonian
		eig_result res = eig_sparse_bandmat(T, max_n, EV_SMALLEST_RE);

		// Check results
		for (u32 n = 0; n < res.num_eigenpairs; ++n) {
			asserteq(f64_compare(res.eigenvalues[n], ho_eigenvalue((i32[]){n}, 1), 1e-5), true);
		}
	}

	it ("HO basis hamiltonian -- pob") {
		// Check results
		//for (u32 n = 0; n < res.num_eigenpairs; ++n) {
		//	c64_normalize(&res.eigenvectors[n*res.points_per_eigenvector], res.points_per_eigenvector);
		//	for (u32 i = 0; i < res.points_per_eigenvector; ++i) {
		//		f64 diff = cabs(res.eigenvectors[n*res.points_per_eigenvector + i]) - fabs(exact_pob_func[n][i]);
		//		asserteq(fabs(diff) < 0.05, true);
		//	}
		//}

		hermitian_bandmat T = construct_ho_kinetic_matrix(max_n + 250);

		u32 params[2];
		integration_settings settings = {
			.gk = gk7,
			.abs_error_tol = 1e-10,
			.rel_error_tol = 1e-10,
			.max_evals = 1e4,
			.userdata = params
		};

		for (u32 r = 0; r < T.size; ++r) {
			for (u32 c = r; c < T.size; ++c) {
				params[0] = r;
				params[1] = c;
				integration_result res = quadgk(pob_integrand, -INFINITY, INFINITY, settings);
				asserteq(res.converged, true);

				u32 i = (T.size-1)*(T.size-(c-r)) + r;
				T.base.data[i] += res.integral;
			}
		}

		asserteq(mat_is_valid(T.base), true);

		// Solve eigenvalue problem for hamiltonian
		eig_result res = eig_sparse_bandmat(T, max_n, EV_SMALLEST_RE);

		{
			plot_init(800, 600, "HO - pob");
			f32 pdata[sample_points];
			c64 cdata[sample_points];
			sample_space sp = make_linspace(1, -L/2.0, L/2.0, sample_points);

			for (u32 i = 0; i < sample_points; ++i) {
				pdata[i] = potential_well_smooth(sp.points[i], L);
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = pdata,
					.label = "potential",
					});

			for (u32 i = 0; i < res.num_eigenpairs; ++i) {
				memset(cdata, 0, sizeof(cdata));
				for (u32 j = 0; j < res.points_per_eigenvector; ++j) {
					c64 coeff = res.eigenvectors[i*res.num_eigenpairs + j];

					// Get the j:th basis function
					for (u32 k = 0; k < sp.pointcount; ++k) {
						f64 x = sp.points[k];
						cdata[k] += coeff*ho_eigenfunction((i32[]){j}, &x, 1);
					}
				}

				for (u32 j = 0; j < sp.pointcount; ++j) {
					pdata[j] = 100*cabs(cdata[j])*cabs(cdata[j]);
				}

				push_line_plot(&(plot_push_desc){
						.space = &sp,
						.data = pdata,
						.label = plot_snprintf("numerical %d", i),
						});
			}
			plot_update_until_closed();
			plot_shutdown();
		}
	}
}
