#include <math.h>
#include <sbmf/common/matrix.h>
#include <sbmf/common/eigenproblem.h>
#include <sbmf/basis/harmonic_oscillator.h>
#include <sbmf/groundstate_solver/groundstate_solver.h>
#include <plot/plot.h>

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

c64 guess(f64* v, i32 n) {
	return gaussian(*v, 0, 0.2);
}

f64 pot(f64* v, i32 n, c64 u) {
	return ho_perturbed_potential(v, n, NULL);
}

describe(item) {
	it ("woo") {
		sbmf_init();
		const f64 L = 5.0;
		const i32 N = 128;
		grid space = generate_grid(1,
				(f64[]){-L/2.0},
				(f64[]){+L/2.0},
				(i32[]){N});
		gss_settings settings = {
			.g = space,
			.max_iterations = 1e4,
			.error_tol = 1e-10,
			.dt = 0.01,
		};
		gss_result item_res = item_execute(settings, pot, guess);
		log_info("item:\niterations: %d\nerror: %e", item_res.iterations, item_res.error);
		c64_normalize(item_res.wavefunction, space.total_pointcount);

		gss_result hob_res = hob(settings, pot, guess);
		log_info("hob:\niterations: %d\nerror: %e", hob_res.iterations, hob_res.error);
		c64_normalize(hob_res.wavefunction, space.total_pointcount);

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
				c64 c = cabs(hob_sample(hob_res.wavefunction, 64, sp.points[i]));
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

		sbmf_shutdown();
	}
}
