#include <sbmf/sbmf.h>
#include <sbmf/methods/gp2c.h>
#include <sbmf/methods/best_meanfield.h>
#include <sbmf/math/functions.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/math/find_eigenpairs.h>

#include <stdio.h>
#include <stdlib.h>

//#define PERTURBATION(x) gaussian(x, 0, 0.2)
#define PERTURBATION(x) (-1.5015*sqrt(x*x - 1.5*1.5 + 1.5015*1.5015));

struct op_params {
	f64 g;
	u32 particle_count;
};

void op(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count],
								void* userdata) {
	struct op_params* params = userdata;
    for (u32 i = 0; i < len; ++i) {
        f64 ca = fabs(in_u[i]);
        out[i] = params->g*(params->particle_count-1)*ca*ca;
    }
}

void perturbation(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count],
								void* userdata) {
    assert(component_count == 0);
    for (u32 i = 0; i < len; ++i) {
        out[i] = PERTURBATION(in_x[i]);
    }
}














void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}





void usage(FILE* stream) {
	fprintf(stream, "Usage: ./bestmf_fig1 [options]\n");
	fprintf(stream, "    -g    specify interaction constant\n");
	fprintf(stream, "    -o    specify output file\n");
}

const char* shiftarg(int* argc, char*** argv) {
	const char* str = **argv;
	(*argc)--;
	(*argv)++;
	return str;
}

int main(int argc, char** argv) {
	u32 particle_count = 100;
	f64 g = (-3.5)/(particle_count-1);
	const char* output_file = "bestmf_output";

	/* skip program name */
	shiftarg(&argc, &argv);

	while (argc > 0) {
		const char* flag = shiftarg(&argc, &argv);
		if (strcmp(flag, "-g") == 0) {
			if (argc == 0) {
				fprintf(stderr, "Error: value for flag %s not specified!\n", flag);
				usage(stderr);
				return 1;
			}

			const char* value = shiftarg(&argc, &argv);
			g = atof(value)/(particle_count - 1);
		} else if (strcmp(flag, "-o") == 0) {
			if (argc == 0) {
				fprintf(stderr, "Error: filename not specified for flag %s!\n", flag);
				usage(stderr);
				return 1;
			}

			output_file = shiftarg(&argc, &argv);
		} else if (strcmp(flag, "-h") == 0) {
			usage(stdout);
			return 0;
		} else {
			fprintf(stderr, "Error: unknown flag %s!\n", flag);
			return 1;
		}
	}

	printf("using\n");
	printf("  particle count: %u\n", particle_count);
	printf("  interaction strength: %lf\n", g);

	sbmf_set_log_callback(log_callback);
	sbmf_init();

	struct gp2c_settings settings = {
        .num_basis_functions = 16,
        .max_iterations = 1e8,
        .error_tol = 1e-9,
        .ho_potential_perturbation = perturbation,
        .gk = gk15,
		.basis = ho_basis,
		.zero_threshold = 1e-10,
    };

	struct op_params params = {
		.g = g,
		.particle_count  = particle_count
	};

    struct gp2c_component components[1] = {
        [0] = {
            .op = op,
			.userdata = &params,
        },
    };

	struct gp2c_result res = gp2c(settings, 1, components);

	u32 bmf_occupation[2];
	f64 bmf_coupling[2*2] = {
		g,g,
		g,g
	};
	f64 bmf_state_coeff[2*res.coeff_count];

	const u32 states_to_include = 2;

	struct eigen_result_real eres = find_eigenpairs_sparse_real(res.hamiltonian[0], states_to_include, EV_SMALLEST_MAG);

	for (u32 i = 0; i < states_to_include; ++i) {
		f64_normalize(&eres.eigenvectors[i*res.coeff_count], &eres.eigenvectors[i*res.coeff_count], res.coeff_count);
		memcpy(&bmf_state_coeff[i*res.coeff_count], &eres.eigenvectors[i*res.coeff_count], res.coeff_count*sizeof(f64));
	}

	FILE* datafile = fopen(output_file, "w");
	for (u32 na0 = 0; na0 <= particle_count; na0 += 5) {
		u32 na1 = particle_count-na0;

		bmf_occupation[0] = na0;
		bmf_occupation[1] = na1;

		f64 E = best_meanfield_energy_new(
				2,
				bmf_occupation,
				bmf_coupling,
				res.coeff_count,
				bmf_state_coeff);

		fprintf(datafile, "%lf\t%lf\n",
				(f64)na0/(f64)particle_count,
				E/(f64)(particle_count)
			   );
	}
	fclose(datafile);


	sbmf_shutdown();
}
