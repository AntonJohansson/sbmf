#include <sbmf/sbmf.h>
#include <stdio.h>

#define PERTURBATION(x) 2*gaussian(x, 0, 0.2)

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


int main() {
	sbmf_set_log_callback(log_callback);

	struct nlse_guess* guesses = NULL;

	u32 occupations[] = {2};

	struct nlse_settings settings = {
		.spatial_pot_perturbation = perturbation,
		.max_iterations = 1000,
		.max_integration_evals = 1e5,
		.error_tol = 1e-9,

        .num_basis_funcs = 16,
		.basis = ho_basis,

		.zero_threshold = 1e-10,
		.gk=gk15
    };

	const u32 component_count = 1;

	for (f64 f = -10; f < 5; f += 0.25) {
		f64 g0[] = {f};
		sbmf_init();
		struct nlse_result res = grosspitaevskii(settings, component_count, occupations, guesses, g0);

		f64 Egp = 0;
		if (res.converged)
			Egp = full_energy(settings, res.coeff_count, component_count, res.coeff, occupations, g0);

		{
			FILE* fd = fopen("out", "a");
			fprintf(fd, "%lf\t%lf\n",
					f,
					Egp
					);
			fclose(fd);
		}

		sbmf_shutdown();
	}

}
