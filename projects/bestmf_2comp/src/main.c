#include <sbmf/sbmf.h>

#include <plot/plot.h>

#include <stdio.h>

#define NA 4
#define NB 4

//#define GAA (1.0/3.0)
#define GAA (-8.0/((f64)NA-1))
#define GAB (+1.0/((f64)NB))
#define GBA (+1.0/((f64)NA))
#define GBB (-8.0/((f64)NB-1))

#define USE_GAUSSIAN_GUESS 0

#define PERTURBATION(x) (cos(5*x))
//#define PERTURBATION(x) 0.0
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
		out[i] = gaussian(in[i] + 1.0, 0.0, 0.2);
}
void gaussian1(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i)
		out[i] = gaussian(in[i] - 2.0, 0.0, 0.1);
}
void gaussian2(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i)
		out[i] = gaussian(in[i] - 1.0, 0.0, 0.2);
}
void gaussian3(f64* out, f64* in, u32 len, void* data) {
	for (u32 i = 0; i < len; ++i)
		out[i] = gaussian(in[i] + 2.0, 0.0, 0.1);
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
		[2] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = gaussian2,
		},
		[3] = {
			.type = SPATIAL_GUESS,
			.data.spatial_guess = gaussian3,
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

	struct nlse_settings settings = {
        .spatial_pot_perturbation = perturbation,
		.max_iterations = 1e5,
		.error_tol = 1e-9,

        .num_basis_funcs = 16,
		.basis = ho_basis,

		.zero_threshold = 1e-13,
		.measure_every = 0,
		.gk=gk15
    };

	const u32 component_count = 2;

	struct bestmf_2comp_result res = best_meanfield_2comp(settings, 4, g0, guesses);
	printf("energy: %lf\n", res.energy);
	printf("n1: %lf\n", res.n1);
	printf("n2: %lf\n", res.n2);
	printf("n3: %lf\n", res.n3);
	printf("n4: %lf\n", res.n4);
	//printf("\nfull energy: %lf\n",
	//		full_energy(settings, res.coeff_count, component_count, res.coeff, occupations, g0));

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

		plot_update_until_closed();
		plot_shutdown();
	}
#endif

	sbmf_shutdown();
}
