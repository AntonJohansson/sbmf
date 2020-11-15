#include <sbmf/sbmf.h>
#include <sbmf/methods/gp2c.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/math/functions.h>
#include <sbmf/methods/best_meanfield.h>
#include <sbmf/math/find_eigenpairs.h>

#include <plot/plot.h>
#define CIMGUI_DEFINE_ENUMS_AND_STRUCTS
#include <cimgui.h>

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#if 1
#define PERTURBATION(x) \
	gaussian(x, 0.0, 0.2)
#else
#define PERTURBATION(x) 0
#endif


#define MAX_COMP_COUNT 16
static f64 gs[MAX_COMP_COUNT*MAX_COMP_COUNT] = {0};
static f64 xoffset[MAX_COMP_COUNT] = {0};

static const char* guess_items[] = {
	"default",
	"gaussian",
};

struct comp_info;

struct comp_info {
	bool being_edited;
	struct comp_info* next;
	u32 cur_guess;

	u32 particle_count;

	/* per particle energies */
	f64 gp_chem_pot;
	f64 gp_full_energy;

	struct eigen_result_real excited_states;
};

static u32 comp_count = 0;
static struct comp_info* comp_info_head = NULL;

static i32 ui_num_basis_funcs = 16;
static f64 ui_error_tol = 1e-9;
static i32 ui_max_iterations = 1e8;
static i32 ui_measure_every = 0;
static f64 ui_zero_threshold = 1e-10;
static i32 ui_excited_states = 2;

/*
 * SBMF
 */

static sample_space sp = {0};
#define L 10.0
#define N 128
static f32 data[N];

static bool gp2c_doing_calc = false;
static bool gp2c_do_plot = false;
static bool gp2c_have_groundstate = false;
static struct gp2c_result gp2c_res;
static struct gp2c_settings gp2c_settings;

static pthread_t thread;
static pthread_mutex_t mutex;





static bool has_excited_states = false;
static i32 ui_states_to_include = 5;




void log_callback(enum sbmf_log_level log_level, const char* msg) {
	switch (log_level) {
		case SBMF_LOG_LEVEL_INFO: 		printf("[info] "); break;
		case SBMF_LOG_LEVEL_WARNING: 	printf("[warning] "); break;
		case SBMF_LOG_LEVEL_ERROR: 		printf("[error] "); break;
		case SBMF_LOG_LEVEL_PANIC: 		printf("[panic] "); break;
		default: printf("[other] "); break;
	};
	printf("%s\n", msg);
}

void debug_callback(struct gp2c_settings s, struct gp2c_result r) {
	pthread_mutex_lock(&mutex);
	gp2c_do_plot = true;
	gp2c_settings = s;
	gp2c_res = r;
	pthread_mutex_unlock(&mutex);
}

/*
 * SBMF guess
 */

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

void guess_callback(f64* out, u32 len, u32 component) {
    static struct integration_settings settings = {
        .abs_error_tol = 1e-10,
        .rel_error_tol = 1e-10,
        .max_evals = 1e5,
    };
    struct guess_integrand_params params = {
        .xoffset = xoffset[component],
    };
    settings.userdata = &params;
    for (u32 i = 0; i < len; ++i) {
        params.n = i;
        struct integration_result res = quadgk_vec(guess_integrand, -INFINITY, INFINITY, settings);
        out[i] = res.integral;
    }
    f64_normalize(out, out, len);
}


/*
 * gp2c
 */

void operator(const u32 len, f64 out[static len],
			  f64 in[static len], const u32 component_count,
			  f64 wf[static len*component_count],
			  void* userdata) {
	for (u32 i = 0; i < len; ++i) {
		out[i] = PERTURBATION(in[i]);
		for (u32 j = 0; j < component_count; ++j) {
			f64* g = userdata;
			f64 c = fabs(wf[j*len + i]);
			out[i] += g[j]*c*c;
		}
	}
}

void* find_groundstate_thread(void* params) {
	pthread_mutex_lock(&mutex);
	gp2c_doing_calc = true;
	gp2c_do_plot = false;
	gp2c_have_groundstate = false;
	/* Copy everything over to the settings */
    struct gp2c_settings settings = {
        .num_basis_functions = ui_num_basis_funcs,
        .max_iterations = ui_max_iterations,
        .error_tol = ui_error_tol,
        //.ho_potential_perturbation = perturbation,
        .gk = gk15,
        .basis = ho_basis,
        .dbgcallback = debug_callback,
        .measure_every = ui_measure_every,
		.zero_threshold = ui_zero_threshold,
    };
	/* setup components */
	struct gp2c_component comps[comp_count];
	{
		struct comp_info* p = comp_info_head;
		u32 iter = 0;
		while (p) {
			//for (u32 i = 0; i < comp_count; ++i) {
			//	gs_w_particle_count[iter*MAX_COMP_COUNT] =
			//		gs[iter*MAX_COMP_COUNT] *
			//		((iter == i) ? (p->particle_count-1) : p->particle_count);
			//}

			comps[iter] = (struct gp2c_component) {
				.op = operator,
				.userdata = &gs[iter*MAX_COMP_COUNT],
			};

			switch (p->cur_guess) {
				case 0: break; /* default */
				case 1: comps[iter].guess = guess_callback; break;
				default: break;/* leaving this blank will result in same as 0 above */
			};

			p = p->next;
			iter++;
		}
	}
	pthread_mutex_unlock(&mutex);

    struct gp2c_result res = gp2c(settings, comp_count, comps);

	pthread_mutex_lock(&mutex);
	gp2c_doing_calc = false;
	gp2c_do_plot = true;
	gp2c_have_groundstate = true;
	gp2c_res = res;
	gp2c_settings = settings;
	pthread_mutex_unlock(&mutex);
}

void find_groundstate() {
	pthread_create(&thread, NULL, find_groundstate_thread, NULL);
}

/*
 * UI
 */

void update() {
	static char buf[64];

	igBegin("settings", 0, 0);
		if (igCollapsingHeaderTreeNodeFlags("settings", 0)) {
			if (igInputInt("basis funcs.", &ui_num_basis_funcs, 1, 8, 0)) {
				if (ui_num_basis_funcs < 1)
					ui_num_basis_funcs = 1;
			}
			if (igInputInt("max iters.", &ui_max_iterations, 1, 8, 0)) {
				if (ui_max_iterations < 1)
					ui_max_iterations = 1;
			}
			if (igInputDouble("error tol.", &ui_error_tol, 0.0, 0.0, "%.2e", 0)) {
				if (ui_error_tol < 0.0)
					ui_error_tol = 0.0;
			}
			if (igInputDouble("zero threshold", &ui_zero_threshold, 0.0, 0.0, "%.2e", 0)) {
				if (ui_zero_threshold < 0.0)
					ui_zero_threshold = 0.0;
			}
			if (igInputInt("measure every", &ui_measure_every, 1, 8, 0)) {
				if (ui_measure_every < 0)
					ui_measure_every = 0;
			}
			igSpacing();
			igSeparator();
			igSpacing();
		}

		/* add comp */
		if (comp_count < MAX_COMP_COUNT && igButton("add component", (ImVec2){0,0})) {
			comp_count++;

			struct comp_info** pp = &comp_info_head;
			while (*pp) {
				pp = &(*pp)->next;
			}

			*pp = malloc(sizeof(struct comp_info));
			memset(*pp, 0, sizeof(struct comp_info));
			(*pp)->particle_count = 100;
		}

		/* comp buttons */
		{
			struct comp_info* p = comp_info_head;
			int iter = 0;
			while (p) {
				snprintf(buf, 64, "component %d", iter);
				if (igButton(buf, (ImVec2){0,0})) {
					p->being_edited = true;
				}
				p = p->next;
				iter++;
			}
		}

		/* g matrix */
		if (comp_count > 0) {
			igSpacing();
			igSeparator();
			igSpacing();

			igText("g matrix");
			struct comp_info* p = comp_info_head;
			u32 iter = 0;
			while (p) {
				snprintf(buf, 64, "comp. %u", iter);
				if (igInputScalarN(buf, ImGuiDataType_Double, &gs[iter*MAX_COMP_COUNT], comp_count, 0, 0, "%.1lf", 0)) {
				}
				p = p->next;
				iter++;
			}
		}

		igSpacing();
		igSeparator();
		igSpacing();

		if (gp2c_doing_calc) {
			igPushStyleColorVec4(ImGuiCol_Button, (ImVec4){0.5,0.1,0.1,1});
		} else {
			igPushStyleColorVec4(ImGuiCol_Button, (ImVec4){0.1,0.5,0.1,1});
		}
		if (igButton("find groundstate", (ImVec2){0,0}) && !gp2c_doing_calc) {
			find_groundstate();
		}
		igPopStyleColor(1);

		if (gp2c_have_groundstate) {
			if (igCollapsingHeaderTreeNodeFlags("Excited states", 0)) {
				if (igInputInt("states to include", &ui_states_to_include, 1,8,0)) {
					if (ui_states_to_include < 1)
						ui_states_to_include = 1;
				}
				if (igButton("calc. excited states", (ImVec2){0,0})) {
					struct comp_info* p = comp_info_head;
					u32 index = 0;
					while (p) {
						p->excited_states = find_eigenpairs_sparse_real(gp2c_res.hamiltonian[index], ui_states_to_include, EV_SMALLEST_MAG);
						p = p->next;
						index++;
					}

					/* Plot the result */
					{
						plot_clear();

						for (u32 i = 0; i < N; ++i) {
							f64 x = sp.points[i];
							data[i] = (f32) ho_potential(&x,1,0) + PERTURBATION(x);
						}
						push_line_plot(&(plot_push_desc){
								.space = &sp,
								.data = data,
								.label = "potential",
								});


						f64 sample_in[N];
						for (u32 i = 0; i < N; ++i) {
							sample_in[i] = (f64) sp.points[i];
						}

						struct comp_info* p = comp_info_head;
						u32 iter = 0;
						while (p) {
							for (u32 j = 0; j < p->excited_states.num_eigenpairs; ++j) {
								f64 sample_out[N];
								gp2c_settings.basis.sample(p->excited_states.points_per_eigenvector,
										&p->excited_states.eigenvectors[j*p->excited_states.points_per_eigenvector],
										N, sample_out, sample_in);

								for (u32 i = 0; i < N; ++i) {
									f64 c = fabs(sample_out[i]);
									data[i] = c*c;
								}
								push_line_plot(&(plot_push_desc){
										.space = &sp,
										.data = data,
										.label = plot_snprintf("comp %u state %u", iter, j),
										.offset = p->excited_states.eigenvalues[j],
										});
							}

							p = p->next;
							iter++;
						}

						for (u32 i = 0; i < gp2c_res.component_count; ++i) {

						}
					}

					has_excited_states = true;
				}
			}
			if (igCollapsingHeaderTreeNodeFlags("groundstate calcs.", 0)) {
				if (igButton("calc. gp energy per particle", (ImVec2){0,0})) {
					struct comp_info* p = comp_info_head;
					u32 iter = 0;
					while (p) {
						p->gp_chem_pot = gp2c_res.energy[iter];
						p->gp_full_energy = gp_energy_per_particle(p->particle_count, gs[iter*MAX_COMP_COUNT + iter]/(p->particle_count-1), gp2c_res.coeff_count, &gp2c_res.coeff[iter*gp2c_res.coeff_count]);
						p = p->next;
						iter++;
					}
				}
				if (igButton("find bestmf occupations", (ImVec2){0,0})) {
				}

				struct comp_info* p = comp_info_head;
				u32 iter = 0;
				while (p) {
					igText("[%u] chem pot: %e", iter, p->gp_chem_pot);
					igText("[%u] full eng: %e", iter, p->gp_full_energy);
					p = p->next;
					iter++;
				}
			}
		}

	igEnd();


	/* Comp windows */
	int index = 0;
	for (struct comp_info* p = comp_info_head; p; p = p->next) {
		if (p->being_edited) {
			snprintf(buf, 64, "Editing component %d", index);
			igBegin(buf, 0,0);
				if (igInputInt("particle count", &p->particle_count, 1, 8, 0)) {
					if (p->particle_count < 1)
						p->particle_count = 1;
				}
				igComboStr_arr("guess", &p->cur_guess, guess_items, sizeof(guess_items)/sizeof(guess_items[0]), 0);
				if (p->cur_guess == 1) {
					/* gaussian */
					igInputDouble("xoffset", &xoffset[index], 0.1, 1.0, "%.3lf", 0);

					/* plot the guess */
					{
						plot_clear();

						for (u32 i = 0; i < N; ++i) {
							f64 x = sp.points[i];
							data[i] = (f32) ho_potential(&x,1,0) + PERTURBATION(x);
						}
						push_line_plot(&(plot_push_desc){
								.space = &sp,
								.data = data,
								.label = "potential",
								});

						for (u32 i = 0; i < N; ++i) {
							data[i] = gaussian((f64)sp.points[i], xoffset[index], 0.2);
						}

						push_line_plot(&(plot_push_desc){
								.space = &sp,
								.data = data,
								.label = "gaussian guess",
								});
					}
				}

				if (igButton("close", (ImVec2){0,0})) {
					p->being_edited = false;
				}

				if (igButton("delete", (ImVec2){0,0})) {
					comp_count--;

					int iter = 0;
					struct comp_info* node = comp_info_head;
					struct comp_info* prev = NULL;
					while (iter < index) {
						prev = node;
						node = node->next;
						iter++;
					}

					if (prev)
						prev->next = node->next;
					else
						comp_info_head = node->next;

					free(node);
				}
			igEnd();
		}
		index++;
	}

	/* Plot result data if we have it */
	if (gp2c_do_plot) {
		/* Plot the result */
		plot_clear();

		for (u32 i = 0; i < N; ++i) {
			f64 x = sp.points[i];
			data[i] = (f32) ho_potential(&x,1,0) + PERTURBATION(x);
		}
		push_line_plot(&(plot_push_desc){
				.space = &sp,
				.data = data,
				.label = "potential",
				});


		f64 sample_in[N];
		for (u32 i = 0; i < N; ++i) {
			sample_in[i] = (f64) sp.points[i];
		}

		for (u32 i = 0; i < gp2c_res.component_count; ++i) {
			f64 sample_out[N];
			gp2c_settings.basis.sample(gp2c_res.coeff_count, &gp2c_res.coeff[i*gp2c_res.coeff_count], N, sample_out, sample_in);

			for (u32 i = 0; i < N; ++i) {
				f64 c = fabs(sample_out[i]);
				data[i] = c*c;
			}
			push_line_plot(&(plot_push_desc){
					.space = &sp,
					.data = data,
					.label = plot_snprintf("comp %u", i),
					.offset = gp2c_res.energy[i],
					});
		}

		gp2c_do_plot = false;
	}
}

int main() {
	sbmf_set_log_callback(log_callback);
	sbmf_init();
	plot_init(800, 600, "explore gp");

	pthread_mutex_init(&mutex, NULL);

	sp = make_linspace(1, -L/2.0, L/2.0, N);

	plot_set_update_callback(update);
	plot_update_until_closed();

	pthread_mutex_destroy(&mutex);

	plot_shutdown();
	sbmf_shutdown();
	return 0;
}
