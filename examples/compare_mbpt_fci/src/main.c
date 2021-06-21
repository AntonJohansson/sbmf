#include <sbmf/sbmf.h>

#include <stdio.h>

#define NA 4
#define NB 0

#define GAA (1.0/3.0)
//#define GAA (-2.0/((f64)NA-1))
#define GAB (+1.0/((f64)NB))
#define GBA (+1.0/((f64)NA))
#define GBB (-2.0/((f64)NB-1))

/* Logging function for sbmf to use */
void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}

int main() {
	sbmf_set_log_callback(log_callback);

	/* NULL means use default, that is eigenstates of the basis */
	struct nlse_guess* guesses = NULL;

	f64 g0[] = {
		GAA, GAB,
		GBA, GBB
	};

	/* Particle couht N = 4 */
	i64 N = 4;

	/* Interaction strengths to try */
	f64 g0s[] = {-1.0/6.0, -1.0/3.0, 1.0/3.0, 1.0/6.0};

	/* Basis sizes to try */
	u32 bs[] = {4,8,12,16,24,32,48,64,80};

	struct nlse_settings settings = {
		.max_iterations = 1e5,
		.max_quadgk_iters = 500,
		.abs_error_tol = 1e-14,

		.num_basis_funcs = 16,
		.basis = ho_basis,

		.zero_threshold = 1e-10,
		.gk=gk20
    };

	const u32 component_count = 1;

	for (u32 j = 0; j < sizeof(g0s)/sizeof(g0s[0]); ++j) {
		f64 g0 = g0s[j];

		for (u32 i = 0; i < sizeof(bs)/sizeof(bs[0]); ++i) {
			u32 b = bs[i];
			settings.num_basis_funcs = b;

			sbmf_init();

			/* Solve GP equation and compute energy */
			struct nlse_result res = grosspitaevskii(settings, component_count, &N, guesses, &g0);
			f64 Egp = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, &N, &g0);

			/* Run perturbation theory */
			struct pt_result rs_ptres = rspt_1comp_cuda(&settings, res, 0, g0, N);
			struct pt_result en_ptres = enpt_1comp_cuda(&settings, res, 0, g0, N);

			/* Save output to file */
			{
				char buf[256];
				snprintf(buf, 256, "out_g%.2lf", g0);

				FILE* fd = fopen(buf, "a");
				fprintf(fd, "%u\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
						b,
						Egp,
						rs_ptres.E0+rs_ptres.E1+rs_ptres.E2,
						rs_ptres.E0+rs_ptres.E1+rs_ptres.E2+rs_ptres.E3,
						en_ptres.E0+en_ptres.E1+en_ptres.E2,
						en_ptres.E0+en_ptres.E1+en_ptres.E2+en_ptres.E3
					   );
				fclose(fd);
			}

			sbmf_shutdown();
		}
	}

}
