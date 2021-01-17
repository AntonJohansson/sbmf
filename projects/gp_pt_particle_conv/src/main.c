#include <sbmf/sbmf.h>

#include <plot/plot.h>

#include <stdio.h>

#define NA 4

#define ORDER (1e5)

#define GAA (1.0/ORDER)

#define USE_GAUSSIAN_GUESS 0
#define USE_TF_GUESS 0

static f64 g0[] = { GAA };
static i64 occupations[] = {NA};

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

void tf(f64* out, f64* in, u32 len, void* data) {
	f64 mu = 0.5 * sqrt(3*g0[0]*(occupations[0]-1)/2);
	for (u32 i = 0; i < len; ++i) {
		out[i] = (mu - 0.5*in[i]*in[i])/4.0;
		if (out[i] < 0)
			out[i] = 0;
		out[i] = sqrt(out[i]);
	}
}

void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}

int main() {
	sbmf_set_log_callback(log_callback);

#if USE_GAUSSIAN_GUESS
	struct nlse_guess guesses[] = {
		[0] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = gaussian0,
		},
	};
#elif USE_TF_GUESS
	struct nlse_guess guesses[] = {
		[0] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = tf,
		},
	};
#else
	struct nlse_guess* guesses = NULL;
#endif


	//u32 bs[] = {4,8,12,16,24,32,48,64};
	//u32 os[] = {5,10,50,100,200,300,400,500,600,700,800,900,1000};
	i64 os[] = {ORDER, 2*ORDER, 3*ORDER, 4*ORDER , 5*ORDER, 6*ORDER, 7*ORDER, 8*ORDER, 9*ORDER};
	//i64 os[] = {9*ORDER};
	struct nlse_settings settings = {
        .spatial_pot_perturbation = perturbation,
		.max_iterations = 1000,
		.max_integration_evals = 1e5,
		.error_tol = 1e-10,

        .num_basis_funcs = 64,
		.basis = ho_basis,

		.hamiltonian_mixing = 0.6,

		.zero_threshold = 1e-10,
		.gk=gk15
    };

	const u32 component_count = 1;

	for (u32 i = 0; i < sizeof(os)/sizeof(os[0]); ++i) {
		u32 o = os[i];
		occupations[0] = o;

		sbmf_init();
		struct nlse_result res = grosspitaevskii(settings, component_count, occupations, guesses, g0);
		f64 Egp = full_energy(settings, res.coeff_count, component_count, res.coeff, occupations, g0);

		//struct pt_result ptres = rayleigh_schroedinger_pt_rf(res, 0, g0, occupations);
		struct pt_result ptres = en_pt_rf(res, 0, g0, occupations);
		//printf("E0:          %.15lf\n", ptres.E0);
		//printf("E1:          %.15lf\n", ptres.E1);
		//printf("E2:          %.15lf\n", ptres.E2);
		//printf("E3:          %.15lf\n", ptres.E3);
		//printf("E0+E1:       %.15lf\n", ptres.E0+ptres.E1);
		//printf("E0+E1+E2:    %.15lf\n", ptres.E0+ptres.E1+ptres.E2);
		//printf("E0+E1+E2+E3: %.15lf\n", ptres.E0+ptres.E1+ptres.E2+ptres.E3);

		{
			FILE* fd = fopen("out", "a");
			f64 Ept = ptres.E0+ptres.E1+ptres.E2+ptres.E3;
			fprintf(fd, "%u\t%.10e\t%.10e\t%.10e\n",
					o,
					Egp/(f64)o,
					Ept/(f64)o,
					(Egp - Ept)/(f64)o
					);
			fclose(fd);
		}

		sbmf_shutdown();
	}

}
