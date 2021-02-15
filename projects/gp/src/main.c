#include <sbmf/sbmf.h>

#include <plot/plot.h>

#include <stdio.h>

//#define NA 4
//#define NB 4
// #define GAA (1.0/(3.0))
// #define GAB (-1.0/6.0)
// #define GBA (-1.0/6.0)
// #define GBB (1.0/3.0)

#define NA 4
#define NB 4
#define GAA (1.0/3.0)
#define GAB (-1.0/6.0)
#define GBA (-1.0/6.0)
#define GBB (1.0/3.0)

#define USE_TF_GUESS 0
#define USE_GAUSSIAN_GUESS 0
#define USE_RANDOM_GUESS 0
#define COMPONENT_COUNT 2

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


void debug_callback(struct nlse_settings settings, struct nlse_result res) {
#if 0
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

#if 0
		{
			FILE* fd = fopen("debug_out", "a");
			fprintf(fd, "%u\t%lf\n", res.iterations, res.energy[0]);
			//for (f64 x = -5; x < 5; x += 0.1) {
			//	f64 out;
			//	ho_sample(res.coeff_count, &res.coeff[0], 1, &out, &x);
			//	fprintf(fd, "%lf\t%lf\n", x, out*out);
			//}
			fclose(fd);
		}
#endif
	}

















void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}

void tf(f64* out, f64* in, u32 len, void* data) {
	f64 mu = 0.5 * sqrt(3*GAA*(NA-1)/2);
	for (u32 i = 0; i < len; ++i) { out[i] = (mu - 0.5*in[i]*in[i])/4.0; if (out[i] < 0)
			out[i] = 0;
		out[i] = sqrt(out[i]);
	}
}

void gaussian0(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i) {
		out[i] = gaussian(in[i] + 1.0, 0.0, 0.1);
	}
}
void gaussian1(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i)
		out[i] = gaussian(in[i] - 1.0, 0.0, 0.1);
}


void expnx(f64* out, f64* in, u32 len, void* p) {
	for (u32 i = 0; i < len; ++i) {
		out[i] = exp(-fabs(in[i]));
	}
}


int main() {
	OMEGA = 1.0;
	sbmf_set_log_callback(log_callback);
	sbmf_init();

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
#elif USE_TF_GUESS
	struct nlse_guess guesses[] = {
		[0] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = tf,
		},
		[1] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = tf,
		},
	};
#elif USE_RANDOM_GUESS
	struct nlse_guess guesses[] = {
		[0] = {
			.type = RANDOM_GUESS,
		},
		[1] = {
			.type = RANDOM_GUESS,
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

	struct nlse_settings settings = {
        .spatial_pot_perturbation = perturbation,
		.max_iterations = 1000,
		.max_quadgk_iters = 500,
		.error_tol = 1e-14,

        .num_basis_funcs = 4,
		.basis = ho_basis,

		.zero_threshold = 1e-10,
		.orbital_mixing = 0.0,
		.hamiltonian_mixing = 0.5,
		.mix_until_iteration = 0,

		.measure_every = 0,
		.debug_callback = debug_callback,
		.gk=gk20
    };

	const u32 component_count = COMPONENT_COUNT;

	struct nlse_result res = grosspitaevskii(settings, component_count, occupations, guesses, g0);
	f64 Efull = grosspitaevskii_energy(settings, res.coeff_count, component_count, res.coeff, occupations, g0);
	printf("\nfull energy: %lf\n", Efull);
	printf("\nfull energy per particle: %lf\n", Efull/((f64)NA+(f64)NB));


	nlse_write_to_binary_file("outbin", res);

#if 1
	{
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


		//{
		//	const u32 S = 5;
		//	//struct eigen_result_real states = find_eigenpairs_sparse_real(res.hamiltonian[0], S, EV_SMALLEST);
		//	struct eigen_result_real states = find_eigenpairs_full_real(res.hamiltonian[0]);
		//	for (u32 j = 0; j < S; ++j) {
		//		f64_normalize(&states.eigenvectors[j*res.coeff_count], &states.eigenvectors[j*res.coeff_count], res.coeff_count);

		//		f64 sample_out[N];
		//		ho_sample(res.coeff_count, &states.eigenvectors[j*res.coeff_count], N, sample_out, sample_in);

		//		f32 data[N];
		//		for (u32 k = 0; k < N; ++k) {
		//			data[k] = fabs(sample_out[k])*fabs(sample_out[k]);
		//		}
		//		push_line_plot(&(plot_push_desc){
		//				.space = &sp,
		//				.data = data,
		//				.label = plot_snprintf("states: %u", j),
		//				.offset = states.eigenvalues[j],
		//				});
		//	}
		//}

		{
			FILE* fd = fopen("output/outplotdata", "w");
			f64 sample[N*res.component_count];
			for (u32 i = 0; i < res.component_count; ++i) {
				f64 sample_out[N];
				ho_sample(res.coeff_count, &res.coeff[i*res.coeff_count], N, &sample[i*N], sample_in);
			}

			for (u32 i = 0; i < N; ++i) {
				fprintf(fd, "%lf\t", sample_in[i]);
				for (u32 j = 0; j < res.component_count; ++j) {
					f64 c = fabs(sample[j*N + i]);
					fprintf(fd, "%lf\t", c*c);
				}
				fprintf(fd, "\n");
			}

			fclose(fd);
		}

		plot_update_until_closed();
		plot_shutdown();
	}
#endif

#if 1
	{
		struct pt_result ptres;
		if (component_count == 2)
			ptres = rayleigh_schroedinger_pt_rf_2comp(settings, res, g0, occupations);
		else
			ptres = rayleigh_schroedinger_pt_rf(settings, res, 0, g0, occupations);
		printf("E0:          %.15e\n", ptres.E0);
		printf("E1:          %.15e\n", ptres.E1);
		printf("E2:          %.15e\n", ptres.E2);
		printf("E3:          %.15e\n", ptres.E3);
		printf("E0+E1:       %.15e\n", ptres.E0+ptres.E1);
		printf("E0+E1+E2:    %.15e\n", ptres.E0+ptres.E1+ptres.E2);
		printf("E0+E1+E2+E3: %.15e\n", ptres.E0+ptres.E1+ptres.E2+ptres.E3);

		printf("diff: %.15lf\n", Efull - (ptres.E0+ptres.E1+ptres.E2+ptres.E3));
	}
#endif

#if 1
	{
		struct pt_result ptres;
		if (component_count == 2)
			ptres = en_pt_2comp(settings, res, g0, occupations);
		else
			ptres = en_pt_rf(settings, res, 0, g0, occupations);
		printf("E0:          %.15e\n", ptres.E0);
		printf("E1:          %.15e\n", ptres.E1);
		printf("E2:          %.15e\n", ptres.E2);
		printf("E3:          %.15e\n", ptres.E3);
		printf("E0+E1:       %.15e\n", ptres.E0+ptres.E1);
		printf("E0+E1+E2:    %.15e\n", ptres.E0+ptres.E1+ptres.E2);
		printf("E0+E1+E2+E3: %.15e\n", ptres.E0+ptres.E1+ptres.E2+ptres.E3);

		printf("diff: %.15lf\n", Efull - (ptres.E0+ptres.E1+ptres.E2+ptres.E3));
	}
#endif

	sbmf_shutdown();
}
