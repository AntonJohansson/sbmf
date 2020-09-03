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
#if 0
		u32 max_n = 128;
		u32 max_repeats = 1;
		u32 max_x = 32;

		f64 x[max_x];
		for (u32 i = 0; i < max_x; ++i) {
			x[i] = i;
		}

		PROFILE_BEGIN("cached precomp");
		hob_precompute_coeffs(CACHE_SIZE);
		PROFILE_END("cached precomp");

		log_info("Comparing sum vs cached sum hob");

		/* Cached sum */
		{
			f64 sum = 0.0;
			for (u32 n = 0; n < max_n; ++n) {
				for (u32 j = 0; j < max_x; ++j) {
					PROFILE_BEGIN("cached sum");
					f64 val = ho_eigenfunction_sumstuff_cached(n,x[j]);
					sum += val;
					PROFILE_END("cached sum");
				}
			}
			log_info("cached sum got: %lf\n", sum);
		}

		/* Cached sum vecx */
		{
			f64 sum = 0.0;
			f64 out[max_x];
			for (u32 n = 0; n < max_n; ++n) {
				PROFILE_BEGIN("cached sum vecx");
				memset(out, 0, sizeof(out));
				ho_eigenfunction_sumstuff_cached_xvec(n,x,out,max_x);
				PROFILE_END("cached sum vecx");
				for (u32 j = 0; j < max_x; ++j) {
					sum += out[j];
				}
			}
			log_info("cached sum vecx got: %lf\n", sum);
		}

		/* sum */
		{
			f64 sum = 0.0;
			for (u32 n = 0; n < max_n; ++n) {
				for (u32 j = 0; j < max_x; ++j) {
					PROFILE_BEGIN("sum");
					sum += ho_eigenfunction_sumstuff(n,x[j]);
					PROFILE_END("sum");
				}
			}
			log_info("sum got: %lf\n", sum);
		}

		/* reccurence */
		{
			f64 sum = 0.0;
			for (u32 n = 0; n < max_n; ++n) {
				for (u32 j = 0; j < max_x; ++j) {
					PROFILE_BEGIN("reccurence");
					sum += ho_eigenfunction((i32[]){n},&x[j],1);
					PROFILE_END("reccurence");
				}
			}
			log_info("reccurence got: %lf\n", sum);
		}

		profile_print_results();
		profile_clear();
#endif

		/* Testing for same value */
		{

			f64 x = 3.0;
			u32 n = 127;
			f64 s1 = ho_eigenfunction((i32[]){n},&x,1);
			f64 s2 = ho_eigenfunction_sumstuff(n,x);
			if (!f64_compare(s1,s2, 0.001)) {
				log_info("failure for n:%u,x=%lf -- s1:%lf, s2:%lf", n,x,s1,s2);
			}

			//for (u32 n = 0; n < max_n; ++n) {
			//	for (u32 j = 0; j < max_x; ++j) {
			//		f64 s1 = ho_eigenfunction((i32[]){n},&x[j],1);
			//		f64 s2 = ho_eigenfunction_sumstuff(n,x[j]);
			//		if (!f64_compare(s1,s2, 0.001)) {
			//			log_info("failure for n:%u,x[%u]=%lf -- s1:%lf, s2:%lf", n,j,x[j],s1,s2);
			//		}
			//	}
			//}
		}
	}
}
