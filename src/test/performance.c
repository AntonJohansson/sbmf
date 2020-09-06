#define SNOW_ENABLED
#include "snow.h"

#include <sbmf/sbmf.h>
#include <sbmf/methods/quadgk.h>
#include <sbmf/debug/profile.h>
#include <sbmf/debug/log.h>

/*
Testing: sin, 0 -> pi: [info] sum: -0.000000 -- iters: 51
[info] Timing results:
[info]                                    entry |  avg. (ns) |   min (ns) |   max (ns) |        tot (ms)
[info] -----------------------------------------+------------+------------+------------+----------------
[info]                                      sin |       1837 |       1818 |       1885 |   32.46966
[info]                    quadgk -- hadapt iter |       1228 |       1214 |       1267 |   21.70652
[info]                      quadgk -- hadapt pq |         53 |         51 |         57 |    0.93954
[info]                    quadgk -- hadapt eval |        988 |        977 |       1029 |   17.47274
✓ Success: sin, 0 -> pi (33.45ms)
quadgk_perf: Passed 1/1 tests. (33.58ms)
*/

void func(f64* out, f64* xvals, u32 len, void* p) {
	for (u32 i = 0; i < len; ++i) {
		out[i] = pow(sin(xvals[i]) + cos(xvals[i]), 3);
	}
}

static void check_quadgk_converge(integration_result res, f64 expected) {
	bool correct_ans = f64_compare(res.integral, expected, 1e-9);
	if (!res.converged || !correct_ans) {
		printf("Integral failed to converge or got wrong answer:\n\tconverged: %d\n\tintegral: %lf\n\terror: %lf\n\texpected: %lf\n", res.converged, res.integral, res.error, expected);
	}

	asserteq(correct_ans && res.converged, true);
}

describe (quadgk_perf) {
	before_each() { sbmf_init(); }
	after_each() { sbmf_shutdown(); }

	integration_settings settings = {
		.gk = gk7,
		.abs_error_tol = 1e-10,
		.rel_error_tol = 1e-10,
		.max_evals = 500,
	};

	it ("sin, 0 -> pi") {
		integration_result res;
		f64 sum = 0.0;
		for (u32 i = 0; i < 10000; ++i) {
			PROFILE_BEGIN("sin");
			res = quadgk(func, -2*M_PI, 2*M_PI, settings);
			PROFILE_END("sin");
			sum += res.integral;
		}
		log_info("sum: %lf -- iters: %d", sum, res.performed_evals);

		profile_print_results();
	}
}

snow_main();
