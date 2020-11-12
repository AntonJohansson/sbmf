#include <sbmf/sbmf.h>
#include <sbmf/methods/gp2c.h>
#include <sbmf/methods/best_meanfield.h>
#include <sbmf/math/functions.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/math/find_eigenpairs.h>

#include <plot/plot.h>

#include <stdio.h>

#define PARTICLE_COUNT 100
#define GAA (-4.0/((f64)PARTICLE_COUNT-1))
#define GAB (+2.0/((f64)PARTICLE_COUNT))
#define GBB (-4.0/((f64)PARTICLE_COUNT-1))
#define PERTURBATION(x) gaussian(x, 0, 0.2)
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

void guess_a(f64* out, u32 len) {
    static struct integration_settings settings = {
        .abs_error_tol = 1e-10,
        .rel_error_tol = 1e-10,
        .max_evals = 1e5,
    };
    struct guess_integrand_params params = {
        .xoffset = 0.25,
    };
    settings.userdata = &params;
    for (u32 i = 0; i < len; ++i) {
        params.n = i;
        struct integration_result res = quadgk_vec(guess_integrand, -INFINITY, INFINITY, settings);
        out[i] = res.integral;
    }
	f64_normalize(out, out, len);
}

void guess_b(f64* out, u32 len) {
    static struct integration_settings settings = {
        .abs_error_tol = 1e-10,
        .rel_error_tol = 1e-10,
        .max_evals = 1e5,
    };
    struct guess_integrand_params params = {
        .xoffset = -0.25,
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
                                f64 in_u[static len*component_count]) {
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

	for (u32 i = 0; i < res.component_count; ++i) {
		printf("[%u] gp single particle energy: %lf\n", i, res.energy[i]);

		f64 E = gp_energy_per_particle(PARTICLE_COUNT, GAA, res.coeff_count, &res.coeff[i*res.coeff_count]);
		printf("[%u] full energy per particle: %lf\n", i, E);
	}

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

		u32 bmf_occupation[4];
		f64 bmf_coupling[4*4] = {
			GAA,GAA,GAB,GAB,
			GAA,GAA,GAB,GAB,
			GBB,GBB,GAB,GAB,
			GBB,GBB,GAB,GAB
		};
		//f64 bmf_coupling[2*2] = {
		//	GAA,GAA,
		//	GAA,GAA,
		//};
		f64 bmf_state_coeff[4*res.coeff_count];

		const u32 states_to_include = 5;
		for (u32 i = 0; i < res.component_count; ++i) {
			struct eigen_result_real eres = find_eigenpairs_sparse_real(res.hamiltonian[i], states_to_include, EV_SMALLEST_MAG);

			for (u32 j = 0; j < states_to_include; ++j) {
				f64_normalize(&eres.eigenvectors[j*res.coeff_count], &eres.eigenvectors[j*res.coeff_count], res.coeff_count);

				f64 sample_out[N];
				ho_sample(res.coeff_count, &eres.eigenvectors[j*res.coeff_count], N, sample_out, sample_in);

				f32 data[N];
				for (u32 k = 0; k < N; ++k) {
					data[k] = fabs(sample_out[k])*fabs(sample_out[k]);
				}
				push_line_plot(&(plot_push_desc){
						.space = &sp,
						.data = data,
						.label = plot_snprintf("comp: %u -- %u", i,j),
						.offset = eres.eigenvalues[j],
						});
			}

			memcpy(&bmf_state_coeff[(2*i+0)*res.coeff_count], &eres.eigenvectors[0*res.coeff_count], res.coeff_count*sizeof(f64));
			memcpy(&bmf_state_coeff[(2*i+1)*res.coeff_count], &eres.eigenvectors[1*res.coeff_count], res.coeff_count*sizeof(f64));

		}

		//f64 best_E = INFINITY;
		//u32 best_na0 = 0;
		//u32 best_nb0 = 0;
		//FILE* datafile = fopen("output/bestmf_data", "w");
		//for (u32 na0 = 0; na0 <= PARTICLE_COUNT; na0 += 5) {
		//	for (u32 nb0 = 0; nb0 <= PARTICLE_COUNT; nb0 += 5) {
		//		//u32 nb0 = 0;
		//		u32 na1 = PARTICLE_COUNT-na0;
		//		u32 nb1 = PARTICLE_COUNT-nb0;

		//		bmf_occupation[0] = na0;
		//		bmf_occupation[1] = na1;
		//		bmf_occupation[2] = nb0;
		//		bmf_occupation[3] = nb1;

		//		f64 E = best_meanfield_energy_new(
		//				4,
		//				bmf_occupation,
		//				bmf_coupling,
		//				res.coeff_count,
		//				bmf_state_coeff);

		//		//f64 E2 = best_meanfield_energy( res.coeff_count,
		//		//							   &bmf_state_coeff[0],
		//		//							   &bmf_state_coeff[res.coeff_count],
		//		//							   na0, na1,
		//		//							   GAA);

		//		fprintf(datafile, "%lf\t%lf\t%lf\n",
		//				(f64)na0/(f64)PARTICLE_COUNT,
		//				(f64)nb0/(f64)PARTICLE_COUNT,
		//				E/(f64)(2.0*PARTICLE_COUNT)
		//				//,E2/(f64)(PARTICLE_COUNT)
		//				);

		//		if (E < best_E) {
		//			best_E = E;
		//			best_na0 = na0;
		//			best_nb0 = nb0;
		//		}
		//	}
		//}
		//fclose(datafile);

		//printf("Best energy: %lf for [%u,%u]\n", best_E/(2.0*PARTICLE_COUNT), best_na0, best_nb0);

		plot_update_until_closed();
		plot_shutdown();
	}



	sbmf_shutdown();
}
