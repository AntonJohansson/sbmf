#include <sbmf/sbmf.h>

#include <plot/plot.h>

#include <stdio.h>

#define NA 4
#define NB 0

#define GAA (1.0/3.0)
//#define GAA (-2.0/((f64)NA-1))
#define GAB (+1.0/((f64)NB))
#define GBA (+1.0/((f64)NA))
#define GBB (-2.0/((f64)NB-1))

#define USE_GAUSSIAN_GUESS 0

//#define PERTURBATION(x) 2*gaussian(x, 0, 0.2)
#define PERTURBATION(x) 0.0
//#define PERTURBATION(x) (-1.5015*sqrt(x*x - 1.5*1.5 + 1.5015*1.5015));

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
		out[i] = gaussian(in[i] + 1.0, 0.0, 0.1);
}
void gaussian1(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i)
		out[i] = gaussian(in[i] - 1.0, 0.0, 0.1);
}




int main() {
	sbmf_set_log_callback(log_callback);

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

	i64 occupations[] = {NA,NB};

	u32 bs[] = {4,8,12,16,24,32,48,64};
	struct nlse_settings settings = {
        .spatial_pot_perturbation = perturbation,
		.max_iterations = 1e5,
		.max_integration_evals = 1e5,
		.error_tol = 1e-9,

        .num_basis_funcs = 16,
		.basis = ho_basis,

		.zero_threshold = 1e-10,
		.gk=gk15
    };

	const u32 component_count = 1;

	for (u32 i = 0; i < sizeof(bs)/sizeof(bs[0]); ++i) {
		u32 b = bs[i];
		settings.num_basis_funcs = b;

		sbmf_init();
		struct nlse_result res = grosspitaevskii(settings, component_count, occupations, guesses, g0);
		f64 Egp = full_energy(settings, res.coeff_count, component_count, res.coeff, occupations, g0);

		//struct pt_result ptres = rayleigh_schroedinger_pt_rf(res, 0, g0, occupations);
		struct pt_result ptres = en_pt_rf(res, 0, g0, occupations);
		printf("E0:          %.15lf\n", ptres.E0);
		printf("E1:          %.15lf\n", ptres.E1);
		printf("E2:          %.15lf\n", ptres.E2);
		printf("E3:          %.15lf\n", ptres.E3);
		printf("E0+E1:       %.15lf\n", ptres.E0+ptres.E1);
		printf("E0+E1+E2:    %.15lf\n", ptres.E0+ptres.E1+ptres.E2);
		printf("E0+E1+E2+E3: %.15lf\n", ptres.E0+ptres.E1+ptres.E2+ptres.E3);

		{
			FILE* fd = fopen("out", "a");
			fprintf(fd, "%u\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
					b,
					Egp, ptres.E0, ptres.E1, ptres.E2, ptres.E3,
					ptres.E0+ptres.E1+ptres.E2+ptres.E3);
			fclose(fd);
		}

		sbmf_shutdown();
	}

}
