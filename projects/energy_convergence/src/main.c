#include <sbmf/sbmf.h>
#include <sbmf/methods/grosspitaevskii.h>
#include <sbmf/methods/best_meanfield.h>
#include <sbmf/math/functions.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/math/find_eigenpairs.h>

#include <plot/plot.h>

#include <stdio.h>
#include <string.h>
#include <time.h>

#define NA 4
#define NB 4

#define GAA (-2.00/((f64)NA-1))
#define GAB (+1.00/((f64)NB))
#define GBA (+1.00/((f64)NA))
#define GBB (-2.00/((f64)NB-1))

#define USE_GAUSSIAN_GUESS 1

#define PERTURBATION(x) 2*gaussian(x, 0, 0.2)
//#define PERTURBATION(x) 0.0
//#define PERTURBATION(x) (-1.5015*sqrt(x*x - 1.5*1.5 + 1.5015*1.5015));

static f64 elapsed_time(struct timespec t0, struct timespec t1) {
	u64 elapsed_ns = (t1.tv_nsec - t0.tv_nsec) + (t1.tv_sec - t0.tv_sec)*(u64)1e9;
	return elapsed_ns/((f64)1e9);
}

static struct timespec current_time() {
	struct timespec t;
	if (clock_gettime(CLOCK_REALTIME, &t) != 0) {
		fprintf(stderr, "Error: Failed to get time!");
	}
	return t;
}

void perturbation(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count],
								void* userdata) {
    for (u32 i = 0; i < len; ++i) {
        out[i] = PERTURBATION(in_x[i]);
    }
}





















void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}




void gaussian0(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i)
		out[i] = gaussian(in[i] + 1.0, 0.0, 0.2);
}
void gaussian1(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i)
		out[i] = gaussian(in[i] - 1.0, 0.0, 0.2);
}




int main() {
	sbmf_set_log_callback(log_callback);

	struct nlse_settings settings = {
		.max_iterations = 1e5,
		.error_tol = 1e-10,
        .spatial_pot_perturbation = perturbation,

		.gk = gk15,

        .num_basis_funcs = 4,
		.basis = ho_basis,

//#if USE_GAUSSIAN_GUESS
//		.guesses = (struct nlse_guess[]) {
//			[0] = {
//				.type = SPATIAL_GUESS,
//				.data.spatial_guess = gaussian0,
//			},
//			[1] = {
//				.type = SPATIAL_GUESS,
//				.data.spatial_guess = gaussian1,
//			},
//		},
//#endif

		.zero_threshold = 1e-10,
    };


#if USE_GAUSSIAN_GUESS
    struct nlse_guess guesses[] = {
        [0] = {
            .type = SPATIAL_GUESS,
            .data.spatial_guess = gaussian0,
        },
        [1] = {
            .type = SPATIAL_GUESS,
            .data.spatial_guess = gaussian1,
        },
    };
#else
    struct nlse_guess* guesses = NULL;
#endif


	f64 g0[] = {
		GAA, GAB,
		GBA, GBB
	};
	u32 occupations[] = {NA,NB};

	u32 comp_count = 2;

	u32 bs[] = {4, 8, 12, 16, 24, 32, 48, 64, /*96,*/ /*128*/};

	for (u32 i = 0; i < sizeof(bs)/sizeof(bs[0]); ++i) {
		FILE* fd = fopen("out", "a");
		sbmf_init();
			u32 b = bs[i];
			printf("b: %u\n", b);
			settings.num_basis_funcs = b;

			struct timespec t0 = current_time();

			struct nlse_result res = grosspitaevskii(settings,
													comp_count,
													occupations,
													guesses,
													g0);

			struct timespec t1 = current_time();

			f64 E = full_energy_naked(settings,
					res.coeff_count, comp_count,
					res.coeff,
					occupations,
					g0);

			fprintf(fd, "%u\t%lf\t%lf\n",
					b,
					E/(occupations[0] + occupations[1]),
					elapsed_time(t0,t1)/res.iterations
					);
		sbmf_shutdown();
		fclose(fd);
	}


}
