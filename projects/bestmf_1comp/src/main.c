#include <sbmf/sbmf.h>

#include <plot/plot.h>

#include <stdio.h>

#define USE_GAUSSIAN_GUESS 1
#define USE_RANDOM_GUESS 0

#define PERTURBATION(x) (-1.5015*sqrt(x*x - 1.5*1.5 + 1.5015*1.5015));
//#define PERTURBATION(x) 2*gaussian(x,0,0.2)
//#define PERTURBATION(x) 0.0

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

	struct nlse_settings settings = {
        .spatial_pot_perturbation = perturbation,
		.max_iterations = 1e5,
		.max_integration_evals = 1e5,
		.error_tol = 1e-9,

        .num_basis_funcs = 32,
		.basis = ho_basis,

		.hamiltonian_mixing = 0.6,
		//.diis_enabled=true,
		//.diis_log_length=8,

		.zero_threshold = 1e-10,
		.measure_every = 0,
		.gk=gk15
    };

	const i64 particle_count = 1000;
	const f64 g0 = (+5.0)/(particle_count-1);

	struct bestmf_result res = best_meanfield(settings, particle_count, g0, guesses);

	//struct nlse_result gp_res_odd = grosspitaevskii(settings, 1, (i64[]){particle_count}, guesses, (f64[]){g0});
	//struct nlse_result gp_res_even = grosspitaevskii(settings, 1, (i64[]){particle_count}, NULL, (f64[]){g0});

	//f64 Egp_odd = full_energy(settings, gp_res_odd.coeff_count, 1, gp_res_odd.coeff, (i64[]){particle_count}, (f64[]){g0});
	//f64 Egp_even = full_energy(settings, gp_res_even.coeff_count, 1, gp_res_even.coeff, (i64[]){particle_count}, (f64[]){g0});
	//printf("gp odd energy: %lf\n", Egp_odd);
	//printf("gp even energy: %lf\n", Egp_even);

	printf("energy: %lf\n", res.energy);
	printf("n1: %lf\n", res.n1);
	printf("n2: %lf\n", res.n2);

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
					//.offset = res.energy[i],
					});
		}

		{
			f64 sample_out1[N];
			ho_sample(res.coeff_count, &res.coeff[0], N, sample_out1, sample_in);
			f64 sample_out2[N];
			ho_sample(res.coeff_count, &res.coeff[res.coeff_count], N, sample_out2, sample_in);

			f32 data[N];
			for (u32 k = 0; k < N; ++k) {
				f64 c1 = fabs(sample_out1[k]);
				f64 c2 = fabs(sample_out2[k]);
				data[k] = (res.n1*c1*c1 + res.n2*c2*c2)/(f64)particle_count;
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = data,
					.label = "bestmf",
					//.offset = res.energy[i],
					});
		}

		/* GP odd */
		//{
		//	f64 sample_out[N];
		//	ho_sample(gp_res_odd.coeff_count, &gp_res_odd.coeff[0], N, sample_out, sample_in);

		//	f32 data[N];
		//	for (u32 k = 0; k < N; ++k) {
		//		data[k] = fabs(sample_out[k])*fabs(sample_out[k]);
		//	}
		//	push_line_plot(&(plot_push_desc){
		//			.space = &sp,
		//			.data = data,
		//			.label = "GPodd",
		//			//.offset = res.energy[i],
		//			});
		//}

		/* GP even */
		//{
		//	f64 sample_out[N];
		//	ho_sample(gp_res_even.coeff_count, &gp_res_even.coeff[0], N, sample_out, sample_in);

		//	f32 data[N];
		//	for (u32 k = 0; k < N; ++k) {
		//		data[k] = fabs(sample_out[k])*fabs(sample_out[k]);
		//	}
		//	push_line_plot(&(plot_push_desc){
		//			.space = &sp,
		//			.data = data,
		//			.label = "GPeven",
		//			//.offset = res.energy[i],
		//			});
		//}

		plot_update_until_closed();
		plot_shutdown();

		/* Save stuff to file */
		{
			f64 out1[N],out2[N],out3[N],out4[N];
			//ho_sample(gp_res_even.coeff_count, &gp_res_even.coeff[0], N, out1, sample_in);
			//ho_sample(gp_res_odd.coeff_count, &gp_res_odd.coeff[0], N, out2, sample_in);
			ho_sample(res.coeff_count, &res.coeff[0], N, out3, sample_in);
			ho_sample(res.coeff_count, &res.coeff[res.coeff_count], N, out4, sample_in);

			FILE* fd = fopen("outplotdata", "w");

			for (u32 k = 0; k < N; ++k) {
				f64 c1 = fabs(out3[k]);
				f64 c2 = fabs(out4[k]);
				f64 ybmf    = (res.n1*c1*c1 + res.n2*c2*c2)/(f64)particle_count;
				//f64 ygpodd  = fabs(out1[k])*fabs(out1[k]);
				//f64 ygpeven = fabs(out2[k])*fabs(out2[k]);

				fprintf(fd, "%lf\t%lf\n", sample_in[k], ybmf);
			}

			fclose(fd);
		}
	}
#endif

	sbmf_shutdown();
}
