#define SNOW_ENABLED
#include "snow.h"

#include <sbmf/sbmf.h>
#include <sbmf/methods/quadgk.h>
#include <sbmf/debug/profile.h>
#include <sbmf/debug/log.h>

f64 func(f64 x, void* p) {
	return pow(sin(x)+cos(x),3);
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
