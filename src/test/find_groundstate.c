#include <math.h>
#include <sbmf/math/matrix.h>
#include <sbmf/math/find_eigenpairs.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/methods/find_groundstate.h>
#include <sbmf/debug/profile.h>
#include <plot/plot.h>

c64 guess(f64* v, i32 n) {
	return gaussian(*v, 0, 0.2);
}

f64 pot(f64* v, i32 n, c64 u) {
	return ho_perturbed_potential(v, n, NULL);
}

describe(item) {
	it ("woo") {
#if 0
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
#endif
	}
}

describe(hob_eigenfunction_perf) {
	before_each(){
		sbmf_init();
	}
	after_each(){
		sbmf_shutdown();
	};

	it ("single x, single n") {
		f64 x = 1.0;

		PROFILE_BEGIN("cached precomp");
		hob_precompute_coeffs(CACHE_SIZE);
		PROFILE_END("cached precomp");

		log_info("Comparing sum vs cached sum hob");
		for (u32 i = 0; i < 500; ++i) {
			for (u32 n = 0; n < 200; ++n) {
				PROFILE_BEGIN("cached sum");
				f64 s1 = ho_eigenfunction_sumstuff_cached(n,x);
				PROFILE_END("cached sum");

				PROFILE_BEGIN("sum");
				f64 s2 = ho_eigenfunction_sumstuff(n,x);
				PROFILE_END("sum");

				PROFILE_BEGIN("reccurence");
				f64 s3 = ho_eigenfunction((i32[]){n},&x,1);
				PROFILE_END("reccurence");

				if (!f64_compare(s1,s2,0.001) || !f64_compare(s1,s3,0.001)) {
					log_error("n: %d -- %lf, %lf, %lf", n,s1,s2,s3);
					assert(false);
				}
			}
		}
	}
	it ("mutliple x, single n") {
	}
	it ("single x, multiple n") {
	}
	it ("multiple x, multiple n") {
	}
}
