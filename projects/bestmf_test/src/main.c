#include <sbmf/sbmf.h>
#include <sbmf/methods/best_meanfield.h>
#include <sbmf/methods/grosspitaevskii.h>
#include <sbmf/math/functions.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/methods/quadgk_vec.h>
#include <plot/plot.h>

#include <stdio.h>

#define PERTURBATION(x) (-1.5015*sqrt(x*x - 1.5*1.5 + 1.5015*1.5015))

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

void gaussian0(f64* out, f64* in, u32 len, void* data) {
    for (u32 i = 0; i < len; ++i)
        out[i] = gaussian(in[i] + 1.0, 0.0, 0.2);
}
void gaussian1(f64* out, f64* in, u32 len, void* data) {
    for (u32 i = 0; i < len; ++i)
        out[i] = gaussian(in[i] - 1.0, 0.0, 0.2);
}

void debug_callback(struct nlse_settings settings, struct nlse_result res) {
#if 1
    const u32 N = 256;
    plot_init(800, 600, "gp2c");

    f32 potdata[N], adata[N], bdata[N];
    sample_space sp = make_linspace(1, -5, 5, N);

    for (u32 i = 0; i < N; ++i) {
        f64 x = sp.points[i];
        potdata[i] = (f32) ho_potential(&x,1,0) + PERTURBATION(x);
    }
    push_line_plot(&(plot_push_desc){
            .space = &sp,
            .data = potdata,
            .label = "potential",
            });

    f64 sample_in[N];
    for (u32 i = 0; i < N; ++i) {
        sample_in[i] = (f64) sp.points[i];
    }

    for (u32 i = 0; i < res.component_count; ++i) {
        f64 sample_out[N];
        ho_sample(res.coeff_count, &res.coeff[i*res.coeff_count], N, sample_out, sample_in);

        f32 data[N];
        for (u32 k = 0; k < N; ++k) {
            data[k] = fabs(sample_out[k])*fabs(sample_out[k]);
        }
        push_line_plot(&(plot_push_desc){
                .space = &sp,
                .data = data,
                .label = plot_snprintf("comp: %u", i),
                .offset = res.energy[i],
                });
    }

    plot_update_until_closed();
    plot_shutdown();
#endif
}

int main(int argc, char** argv) {
	u32 particle_count = 1000;

	sbmf_set_log_callback(log_callback);

	struct nlse_settings settings = {
		.num_basis_funcs = 16,
		.max_iterations = 1e5,
		.error_tol = 1e-9,
		.spatial_pot_perturbation = perturbation,
		.gk = gk15,
		.basis = ho_basis,
		.zero_threshold = 1e-10,
		.debug_callback = debug_callback,
		.measure_every = 3,
	};

	struct nlse_guess gaussian_guesses[2] = {
		[0] = { .type = SPATIAL_GUESS, .data.spatial_guess = gaussian0 },
		[1] = { .type = SPATIAL_GUESS, .data.spatial_guess = gaussian1 },
	};

	struct nlse_guess default_guesses[2] = {
		[0] = { .type = DEFAULT_GUESS },
		[1] = { .type = DEFAULT_GUESS },
	};

#define NUM_GUESSES 2
	struct nlse_guess* guesses_to_try[NUM_GUESSES] = {
		gaussian_guesses,
		default_guesses
	};

	f64 gs[] = {-4.0};//-4.0, -3.5, -2.5, -1.5, -1.0, -0.5, 0.5};
	for (u32 i = 0; i < sizeof(gs)/sizeof(gs[0]); ++i) {
		FILE* fd = fopen("out", "a");
		f64 g = (gs[i])/(particle_count-1);

		/* GP problem */
		settings.measure_every = 0;
		f64 Egp = 0.0;
		{
			sbmf_init();
			struct nlse_result gpres = grosspitaevskii(settings, 1, &particle_count, NULL, &g);
			Egp = full_energy_naked(settings,
								     gpres.coeff_count, 1,
								     gpres.coeff, &particle_count,
								     &g);
			sbmf_shutdown();
		}

		/* Betsmf */
		sbmf_init();


		settings.measure_every = 11;
		fprintf(fd, "%lf", gs[i]);
		for (u32 j = 0; j < NUM_GUESSES; ++j) {
			struct bestmf_result res = best_meanfield(settings, particle_count, g, guesses_to_try[j]);

			f64 Eoffset = res.energy - Egp;

			fprintf(fd, "|\t%lf\t%lf\t%lf\t%lf",
					res.n1/((f64)particle_count),
					res.n2/((f64)particle_count),
					res.energy/(f64)particle_count,
					Eoffset/(f64)particle_count);

#if 1
			{
				plot_init(800,600,"?");

				const u32 N = 512;
				sample_space sp = make_linspace(1, -5, 5, N);

				f64 sample_in[N];
				for (u32 i = 0; i < N; ++i) {
					sample_in[i] = (f64) sp.points[i];
				}

				for (u32 i = 0; i < res.comp_count; ++i) {
					f64 sample_out[N];
					ho_sample(res.coeff_count, &res.coeff[i*res.coeff_count], N, sample_out, sample_in);

					f32 data[N];
					for (u32 k = 0; k < N; ++k) {
						data[k] = fabs(sample_out[k])*fabs(sample_out[k]);
					}
					push_line_plot(&(plot_push_desc){
							.space = &sp,
							.data = data,
							.label = plot_snprintf("comp: %u", i),
							.offset = 0,
							});
				}

				plot_update_until_closed();
				plot_shutdown();
			}
#endif
		}
		fprintf(fd, "\n");




		sbmf_shutdown();

		fclose(fd);
	}
}
