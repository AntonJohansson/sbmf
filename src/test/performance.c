#define SNOW_ENABLED
#include "snow.h"

#include <sbmf/sbmf.h>
#include <sbmf/methods/quadgk.h>
#include <sbmf/methods/quadgk_vec.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/debug/profile.h>
#include <sbmf/debug/log.h>

f64 func(f64 x, void* p) {
	f64 eig1 = ho_eigenfunction((i32[]){1}, &x, 1);
	f64 eig2 = ho_eigenfunction((i32[]){2}, &x, 1);
	f64 pot  = ho_potential(&x, 1, 0);
	return eig1*eig2*pot;
}

void func_vec(f64* out, f64* in, u32 len, void* p) {
	f64 eigs1[len];
	f64 eigs2[len];
	f64 pots[len];
	ho_eigenfunction_vec(1, eigs1, in, len);
	ho_eigenfunction_vec(2, eigs2, in, len);
	ho_potential_vec(pots, in, len);

	for (u32 i = 0; i < len; ++i) {
		out[i] = eigs1[i]*eigs2[i]*pots[i];
	}
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
		const u32 repetitions = 10000;
        integration_result res;

        f64 sum = 0.0;
        for (u32 i = 0; i < repetitions; ++i) {
            PROFILE_BEGIN("quadgk");
            res = quadgk(func, -INFINITY, INFINITY, settings);
            PROFILE_END("quadgk");
			sum += res.integral;
		}
		log_info("sum: %lf -- iters: %d", sum, res.performed_evals);

        f64 sum_vec = 0.0;
        for (u32 i = 0; i < repetitions; ++i) {
            PROFILE_BEGIN("quadgk_vec");
            res = quadgk_vec(func_vec, -INFINITY, INFINITY, settings);
            PROFILE_END("quadgk_vec");
			sum_vec += res.integral;
		}
		log_info("sum_vec: %lf -- iters: %d", sum_vec, res.performed_evals);

		profile_print_results();
	}
}

snow_main();
