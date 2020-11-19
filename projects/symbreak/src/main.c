#include <sbmf/sbmf.h>
#include <sbmf/methods/gp2c.h>
#include <sbmf/methods/best_meanfield.h>
#include <sbmf/math/functions.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/math/find_eigenpairs.h>

#include <plot/plot.h>

#include <stdio.h>

#define PARTICLE_COUNT 1000
#define GAA (+1.0/((f64)PARTICLE_COUNT-1))
#define GAB (+2.0/((f64)PARTICLE_COUNT))
#define GBB (+1.0/((f64)PARTICLE_COUNT-1))
#define PERTURBATION(x) 0.5*gaussian(x, 0, 0.1)
//#define PERTURBATION(x) 0.0
//#define PERTURBATION(x) (-1.5015*sqrt(x*x - 1.5*1.5 + 1.5015*1.5015));

void op_a(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count],
								void* userdata) {
    for (u32 i = 0; i < len; ++i) {
        f64 ca = fabs(in_u[i]);
        f64 cb = fabs(in_u[len + i]);
        out[i] =
            + GAA*(PARTICLE_COUNT-1)*ca*ca
            + GAB*(PARTICLE_COUNT)*cb*cb;
    }
}

void op_b(const u32 len, f64 out[static len],
                                f64 in_x[static len], const u32 component_count,
                                f64 in_u[static len*component_count],
								void* userdata) {
    for (u32 i = 0; i < len; ++i) {
        f64 ca = fabs(in_u[i]);
        f64 cb = fabs(in_u[len + i]);
        out[i] =
            + GBB*(PARTICLE_COUNT-1)*cb*cb
            + GAB*(PARTICLE_COUNT)*ca*ca;
    }
}

struct guess_integrand_params {
    u32 n;
    f64 xoffset;
};

void guess_integrand(f64* out, f64* in, u32 len, void* data) {
    struct guess_integrand_params* params = data;

	f64 eig;
	ho_eigenfunc(params->n, 1, &eig, in);

    for (u32 i = 0; i < len; ++i) {
        out[i] = eig * gaussian(in[i], params->xoffset, 0.2);
    }
}

void guess_a(f64* out, u32 len, u32 comp) {
    static struct integration_settings settings = {
        .abs_error_tol = 1e-10,
        .rel_error_tol = 1e-10,
        .max_evals = 1e5,
    };
    struct guess_integrand_params params = {
        .xoffset = 1,
    };
    settings.userdata = &params;
    for (u32 i = 0; i < len; ++i) {
        params.n = i;
        struct integration_result res = quadgk_vec(guess_integrand, -INFINITY, INFINITY, settings);
        out[i] = res.integral;
    }
	f64_normalize(out, out, len);
}

void guess_b(f64* out, u32 len, u32 comp) {
    static struct integration_settings settings = {
        .abs_error_tol = 1e-10,
        .rel_error_tol = 1e-10,
        .max_evals = 1e5,
    };
    struct guess_integrand_params params = {
        .xoffset = -1,
    };
    settings.userdata = &params;
    for (u32 i = 0; i < len; ++i) {
        params.n = i;
        struct integration_result res = quadgk_vec(guess_integrand, -INFINITY, INFINITY, settings);
        out[i] = res.integral;
    }
	f64_normalize(out, out, len);
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













void debug_callback(struct gp2c_settings settings, struct gp2c_result res) {

  plot_init(800,600,"debug");

  const u32 N = 256;
  f32 adata[N];
  sample_space sp = make_linspace(1, -5, 5, N);

  f64 sample_in[N];
  for (u32 i = 0; i < N; ++i) {
      sample_in[i] = (f64) sp.points[i];
  }

  for (u32 i = 0; i < res.component_count; ++i) {
	  f64 sample_out[N];
	  ho_sample(res.coeff_count, &res.coeff[i*res.coeff_count], N, sample_out, sample_in);

	  for (u32 i = 0; i < N; ++i) {
		  f64 ca = fabs(sample_out[i]);
		  adata[i] = ca*ca;
	  }

	  push_line_plot(&(plot_push_desc){
			  .space = &sp,
			  .data = adata,
			  .label = plot_snprintf("comp %u", i),
			  });
  }


  plot_update_until_closed();
  plot_shutdown();
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









int main() {
	sbmf_set_log_callback(log_callback);
	sbmf_init();

	printf("gab^2 < gaa*gbb | %lf < %lf\n", GAB*GAB, GAA*GBB);

	struct gp2c_settings settings = {
        .num_basis_functions = 16,
        .max_iterations = 1e8,
        .error_tol = 1e-9,
        .ho_potential_perturbation = perturbation,
        .gk = gk15,
		.basis = ho_basis,
		.zero_threshold = 1e-10, /* ? */
		//.dbgcallback = debug_callback,
		//.measure_every = 1,
    };
    struct gp2c_component components[2] = {
        [0] = {
            .op = op_a,
            .guess = guess_a
        },
        [1] = {
            .op = op_b,
            .guess = guess_b
        },
    };

	struct gp2c_result res = gp2c(settings, 2, components);

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

		f64 E = gp_energy_per_particle(PARTICLE_COUNT, GAA, res.coeff_count, &res.coeff[i*res.coeff_count]);
		printf("[%u] full energy per particle: %lf\n", i, E);
	}

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







	f64 pt1_coeffs[res.coeff_count];
	memcpy(pt1_coeffs, res.coeff, res.coeff_count*sizeof(f64));

	printf("Trying perturbation theory for comp 0\n");
	{
		const u32 NA = PARTICLE_COUNT;
		const u32 NB = PARTICLE_COUNT;
		const u32 l = res.coeff_count;

#define PHI(i,j) 	states[i].eigenvectors[j*res.coeff_count]
#define ENERGY(i,j) states[i].eigenvalues[j]

		/* <0|Hmf|0> */
		f64 E0 = NA*ENERGY(0,0) + NB*ENERGY(1,0);

		/* <0|V|0> */

		f64 E1 = 	-0.5*GAA*NA*(NA-1)*V(l, &PHI(0,0), &PHI(0,0), &PHI(0,0), &PHI(0,0))
					-0.5*GBB*NB*(NB-1)*V(l, &PHI(1,0), &PHI(1,0), &PHI(1,0), &PHI(1,0))
					        -GAB*NA*NB*V(l, &PHI(0,0), &PHI(1,0), &PHI(0,0), &PHI(1,0));

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
				E2 += tmp*tmp/(
						NA*ENERGY(0,0)+NB*ENERGY(1,0)
						- ((NA-2)*ENERGY(0,0) + ENERGY(0,i) + ENERGY(0,j))
						- NB*ENERGY(1,0)
						);

				f64 scaling = tmp/(
						NA*ENERGY(0,0)+NB*ENERGY(1,0)
						- ((NA-2)*ENERGY(0,0) + ENERGY(0,i) + ENERGY(0,j))
						- NB*ENERGY(1,0)
						);
				for (u32 k = 0; k < res.coeff_count; ++k) {
					pt1_coeffs[k] += scaling*((&PHI(0,i))[k] + (&PHI(0,j))[k]);
				}
				/* B comp */
				tmp = factor * GBB*sqrt(NB*(NB-1))*V(l, &PHI(1,i), &PHI(1,j), &PHI(1,0), &PHI(1,0));
				E2 += tmp*tmp/(
						NA*ENERGY(0,0)+NB*ENERGY(1,0)
						- NA*ENERGY(0,0)
						- ((NB-2)*ENERGY(1,0) + ENERGY(1,i) + ENERGY(1,j))
						);
			}
		}

		/* Double substitutions (separate components) */
		for (u32 i = 1; i < states_to_include; ++i) {
			for (u32 j = 1; j < states_to_include; ++j) {
				f64 tmp;
				tmp = GAB*sqrt(NA*NB)*V(l, &PHI(0,i), &PHI(1,j), &PHI(0,0), &PHI(1,0));
				E2 += tmp*tmp/(
						NA*ENERGY(0,0) + NB*ENERGY(1,0)
						- ((NA-1)*ENERGY(0,0) + ENERGY(0,i))
						- ((NB-1)*ENERGY(1,0) + ENERGY(1,j))
						);


				f64 scaling = tmp/(
						NA*ENERGY(0,0) + NB*ENERGY(1,0)
						- ((NA-1)*ENERGY(0,0) + ENERGY(0,i))
						- ((NB-1)*ENERGY(1,0) + ENERGY(1,j))
						);
				for (u32 k = 0; k < res.coeff_count; ++k) {
					pt1_coeffs[k] += scaling*((&PHI(0,i))[k]);
				}

			}
		}

		printf("E0: %.2e\t%.2lf\n", E0, E0/(NA+NB));
		printf("E1: %.2e\t%.2lf\n", E1, E1/(NA+NB));
		printf("E2: %.2e\t%.2lf\n", E2, E2/(NA+NB));
		printf("E0+E1+E2 = %.2e\t%.2lf\n", E0+E1+E2, (E0+E1+E2)/(NA+NB));
	}

	f64_normalize(pt1_coeffs, pt1_coeffs, res.coeff_count);

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

		f64 sample_out[N];
		ho_sample(res.coeff_count, pt1_coeffs, N, sample_out, sample_in);

		f32 data[N];
		for (u32 k = 0; k < N; ++k) {
			data[k] = fabs(sample_out[k])*fabs(sample_out[k]);
		}
		push_line_plot(&(plot_push_desc){
				.space = &sp,
				.data = data,
				.label = plot_snprintf("?"),
				});

		plot_update_until_closed();
		plot_shutdown();
	}
#endif

	sbmf_shutdown();
}
