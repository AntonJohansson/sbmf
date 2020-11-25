#include <sbmf/sbmf.h>
#include <sbmf/methods/grosspitaevskii.h>
#include <sbmf/methods/best_meanfield.h>
#include <sbmf/math/functions.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/math/find_eigenpairs.h>

#include <plot/plot.h>

#include <stdio.h>

#define NA 4
#define NB 4

#define GAA (+1/((f64)NA-1))
#define GAB (+2.00/((f64)NB))
#define GBA (+2.00/((f64)NA))
#define GBB (+1/((f64)NB-1))

#define USE_GAUSSIAN_GUESS 1

//#define PERTURBATION(x) gaussian(x, 0, 0.2)
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

	//f32 g[] = {
	//	1, 2,
	//	2, 1
	//};

	//for (u32 i = 0; i < res.component_count; ++i) {
	//	f32 data[N];
	//	memset(data, 0, N*sizeof(f32));

	//	for (u32 j = 0; j < res.component_count; ++j) {
	//		f64 sample_out[N];
	//		ho_sample(res.coeff_count, &res.coeff[j*res.coeff_count], N, sample_out, sample_in);

	//		for (u32 k = 0; k < N; ++k) {
	//			data[k] += g[i*2 + j]*fabs(sample_out[k])*fabs(sample_out[k]);
	//		}
	//	}

	//	push_line_plot(&(plot_push_desc){
	//			.space = &sp,
	//			.data = data,
	//			.label = plot_snprintf("%u ?", i),
	//			.offset = res.energy[i],
	//			});
	//}

	plot_update_until_closed();
	plot_shutdown();
#endif
}

















void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}



struct V_params {
	u32 coeff_count;
	f64* i;
	f64* j;
	f64* k;
	f64* l;
};

void V_integrand(f64* out, f64* in, u32 len, void* data) {
	struct V_params* p = data;

	f64 sample_i[len]; ho_sample(p->coeff_count, p->i, len, sample_i, in);
	f64 sample_j[len]; ho_sample(p->coeff_count, p->j, len, sample_j, in);
	f64 sample_k[len]; ho_sample(p->coeff_count, p->k, len, sample_k, in);
	f64 sample_l[len]; ho_sample(p->coeff_count, p->l, len, sample_l, in);

	for (u32 i = 0; i < len; ++i) {
		out[i] = sample_i[i]*sample_j[i]*sample_k[i]*sample_l[i];
	}
}

f64 V(const u32 coeff_count, f64 i[static coeff_count],
							 f64 j[static coeff_count],
							 f64 k[static coeff_count],
							 f64 l[static coeff_count]) {
	struct V_params p = {
		.coeff_count = coeff_count,
		.i = i,
		.j = j,
		.k = k,
		.l = l
	};

    struct integration_settings settings = {
        .abs_error_tol = 1e-10,
        .rel_error_tol = 1e-10,
        .max_evals = 1e5,
		.userdata = &p,
    };

	struct integration_result res = quadgk_vec(V_integrand, -INFINITY, INFINITY, settings);
	assert(res.converged);

	return res.integral;
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

	printf("gab^2 < gaa*gbb | %lf < %lf\n", GAB*GAB, GAA*GBB);

	struct gp_settings settings = {
        .pot = perturbation,

        .num_basis_funcs = 16,
		.basis = ho_basis,

		.component_count = 2,
		.occupations = (u32[]){NA,NB},
#if USE_GAUSSIAN_GUESS
		.guesses = (struct nlse_guess[]) {
			[0] = {
				.type = SPATIAL_GUESS,
				.data.spatial_guess = gaussian0,
			},
			[1] = {
				.type = SPATIAL_GUESS,
				.data.spatial_guess = gaussian1,
			},
		},
#endif
		.g0 = (f64[]){
			GAA, GAB,
			GBA, GBB
		},
		.zero_threshold = 1e-10,
		.debug_callback = debug_callback,
		.measure_every = 0,
    };

	struct nlse_result res = grosspitaevskii(settings);
	printf("\nfull energy per particle: %lf\n", full_energy(settings, res)/(settings.occupations[0] + settings.occupations[1]));

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

		plot_update_until_closed();
		plot_shutdown();
	}
#endif

	for (u32 i = 0; i < res.component_count; ++i) {
		printf("[%u] gp single particle energy: %lf\n", i, res.energy[i]);
	}
	printf("\n");

	const u32 states_to_include = 5;

	u32 bmf_occupation[4];

	//f64 bmf_coupling[4*4] = {
	//	GAA,GAA,GAB,GAB,
	//	GAA,GAA,GAB,GAB,
	//	GBB,GBB,GAB,GAB,
	//	GBB,GBB,GAB,GAB
	//};

	f64 bmf_coupling[2*2] = {
		GAA,GAA,
		GAA,GAA,
	};

	f64 bmf_state_coeff[4*res.coeff_count];
	f64 bmf_state_energy[4];

	struct eigen_result_real states[res.component_count];

	for (u32 i = 0; i < res.component_count; ++i) {
		states[i] = find_eigenpairs_sparse_real(res.hamiltonian[i], states_to_include, EV_SMALLEST_MAG);

		for (u32 j = 0; j < states_to_include; ++j) {
			f64_normalize(&states[i].eigenvectors[j*res.coeff_count], &states[i].eigenvectors[j*res.coeff_count], res.coeff_count);
			//bmf_state_energy[2*i+j] = eres.eigenvalues[j];
			//memcpy(&bmf_state_coeff[(2*i+j)*res.coeff_count], &eres.eigenvectors[j*res.coeff_count], res.coeff_count*sizeof(f64));
		}
	}

#if 0
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
			for (u32 j = 0; j < states_to_include; ++j) {
				f64 sample_out[N];
				ho_sample(res.coeff_count, &states[i].eigenvectors[j*res.coeff_count], N, sample_out, sample_in);

				f32 data[N];
				for (u32 k = 0; k < N; ++k) {
					data[k] = fabs(sample_out[k])*fabs(sample_out[k]);
				}
				push_line_plot(&(plot_push_desc){
						.space = &sp,
						.data = data,
						.label = plot_snprintf("comp: %u -- %u", i,j),
						.offset = states[i].eigenvalues[j]
						});
			}
		}

		plot_update_until_closed();
		plot_shutdown();
	}
#endif

#if 0
	printf("Trying bestmf for comp 0\n");
	{
		f64 best_E = INFINITY;
		u32 best_na0 = 0;
		FILE* datafile = fopen("output/bestmf_data", "w");
		for (u32 na0 = 0; na0 <= PARTICLE_COUNT; na0 += 10) {
				u32 na1 = PARTICLE_COUNT-na0;

				bmf_occupation[0] = na0;
				bmf_occupation[1] = na1;

				f64 E = best_meanfield_energy_new(
						2,
						bmf_occupation,
						bmf_coupling,
						res.coeff_count,
						bmf_state_coeff);

				fprintf(datafile, "%lf\t%lf\n",
						(f64)na0/(f64)PARTICLE_COUNT,
						E/(f64)(PARTICLE_COUNT)
						);

				if (E < best_E) {
					best_E = E;
					best_na0 = na0;
				}
		}
		fclose(datafile);

		printf("Best energy: %lf for %u\n", best_E/(PARTICLE_COUNT), best_na0);
	}
#endif







	f64 pt1_coeffs_a[res.coeff_count];
	memcpy(pt1_coeffs_a, res.coeff, res.coeff_count*sizeof(f64));

	//f64 pt1_coeffs_b[res.coeff_count];
	//memcpy(pt1_coeffs_b, &res.coeff[res.coeff_count], res.coeff_count*sizeof(f64));

	printf("Trying perturbation theory for comp 0\n");
	{
		const u32 l = res.coeff_count;

#define PHI(i,j) 	states[i].eigenvectors[j*res.coeff_count]
#define ENERGY(i,j) states[i].eigenvalues[j]

		/* <0|Hmf|0> */
		f64 E0 = NA*ENERGY(0,0);// + NB*ENERGY(1,0);

		/* <0|V|0> */

		f64 E1 = 	-0.5*GAA*NA*(NA-1)*V(l, &PHI(0,0), &PHI(0,0), &PHI(0,0), &PHI(0,0))
					;//-0.5*GBB*NB*(NB-1)*V(l, &PHI(1,0), &PHI(1,0), &PHI(1,0), &PHI(1,0))
					 //        -GAB*NA*NB*V(l, &PHI(0,0), &PHI(1,0), &PHI(0,0), &PHI(1,0));

		/* sum |<k|V|0>|^2/(E0-Ek) */
		f64 E2 = 0.0;

		/* Single substitutions are zero!*/

		/* Double substitutions (both excitations within same component) */
		for (u32 i = 1; i < states_to_include; ++i) {
			for (u32 j = i; j < states_to_include; ++j) {
				f64 factor = (i == j) ? 1.0/sqrt(2.0) : 1.0;
				f64 tmp;
				/* A comp */
				tmp = factor * GAA*sqrt(NA*(NA-1))*V(l, &PHI(0,i), &PHI(0,j), &PHI(0,0), &PHI(0,0));
				E2 += tmp*tmp/(2*ENERGY(0,0) - (ENERGY(0,i) + ENERGY(0,j)));

				f64 scaling = tmp/(2*ENERGY(0,0) - (ENERGY(0,i) + ENERGY(0,j)));
				for (u32 k = 0; k < res.coeff_count; ++k) {
					pt1_coeffs_a[k] += scaling*((&PHI(0,i))[k] + (&PHI(0,j))[k]);
					//pt1_coeffs_b[k] += scaling*((&PHI(1,i))[k] + (&PHI(1,j))[k]);
				}
				/* B comp */
				//tmp = factor * GBB*sqrt(NB*(NB-1))*V(l, &PHI(1,i), &PHI(1,j), &PHI(1,0), &PHI(1,0));
				//E2 += tmp*tmp/(
				//		NA*ENERGY(0,0)+NB*ENERGY(1,0)
				//		- NA*ENERGY(0,0)
				//		- ((NB-2)*ENERGY(1,0) + ENERGY(1,i) + ENERGY(1,j))
				//		);
			}
		}

		/* Double substitutions (separate components) */
		//for (u32 i = 1; i < states_to_include; ++i) {
		//	for (u32 j = 1; j < states_to_include; ++j) {
		//		f64 tmp;
		//		tmp = GAB*sqrt(NA*NB)*V(l, &PHI(0,i), &PHI(1,j), &PHI(0,0), &PHI(1,0));
		//		E2 += tmp*tmp/(
		//				NA*ENERGY(0,0) + NB*ENERGY(1,0)
		//				- ((NA-1)*ENERGY(0,0) + ENERGY(0,i))
		//				- ((NB-1)*ENERGY(1,0) + ENERGY(1,j))
		//				);


		//		f64 scaling = tmp/(
		//				NA*ENERGY(0,0) + NB*ENERGY(1,0)
		//				- ((NA-1)*ENERGY(0,0) + ENERGY(0,i))
		//				- ((NB-1)*ENERGY(1,0) + ENERGY(1,j))
		//				);
		//		for (u32 k = 0; k < res.coeff_count; ++k) {
		//			pt1_coeffs_a[k] += scaling*((&PHI(0,i))[k]);
		//			//pt1_coeffs_b[k] += scaling*((&PHI(1,i))[k]);
		//		}

		//	}
		//}

		printf("\nPT energy          PT energy per particle\n");
		printf("E0: %.2e             %lf\n", E0, E0/(NA+NB));
		printf("E1: %.2e             %lf\n", E1, E1/(NA+NB));
		printf("E2: %.2e             %lf\n", E2, E2/(NA+NB));
		printf("E0+E1+E2 = %.2e      %lf\n", E0+E1+E2, (E0+E1+E2)/(NA+NB));
		printf("\n");
	}

	f64_normalize(pt1_coeffs_a, pt1_coeffs_a, res.coeff_count);
	//f64_normalize(pt1_coeffs_b, pt1_coeffs_b, res.coeff_count);

#if 0
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

		f64 sample_a[N];
		ho_sample(res.coeff_count, pt1_coeffs_a, N, sample_a, sample_in);

		//f64 sample_b[N];
		//ho_sample(res.coeff_count, pt1_coeffs_b, N, sample_b, sample_in);

		f32 data_a[N], data_b[N];
		for (u32 k = 0; k < N; ++k) {
			data_a[k] = fabs(sample_a[k])*fabs(sample_a[k]);
			//data_b[k] = fabs(sample_b[k])*fabs(sample_b[k]);
		}
		push_line_plot(&(plot_push_desc){
				.space = &sp,
				.data = data_a,
				.label = plot_snprintf("a ?"),
				});

		//push_line_plot(&(plot_push_desc){
		//		.space = &sp,
		//		.data = data_b,
		//		.label = plot_snprintf("b ?"),
		//		});

		plot_update_until_closed();
		plot_shutdown();
	}
#endif

	sbmf_shutdown();
}
