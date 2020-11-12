#include <sbmf/sbmf.h>
#include <sbmf/methods/gp2c.h>
#include <sbmf/math/harmonic_oscillator.h>
#include <sbmf/math/functions.h>

#include <plot/plot.h>
#define CIMGUI_DEFINE_ENUMS_AND_STRUCTS
#include <cimgui.h>

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

struct component_info;

struct component_info {
	bool being_edited;
	double interaction_strength;
	struct component_info* next;
};

static struct component_info* comp_info_head = NULL;

static u32 ui_num_basis_funcs = 16;
static f64 ui_error_tol = 1e-9;
static u32 ui_max_iterations = 1e8;
static u32 ui_measure_every = 0;

/*
 * SBMF
 */

static sample_space sp = {0};
#define L 10.0
#define N 128
static f32 data[N];

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

void debug_callback(struct gp2c_settings settings, struct gp2c_result res) {
	plot_clear();

	for (u32 i = 0; i < N; ++i) {
		f64 x = sp.points[i];
		data[i] = (f32) ho_potential(&x,1,0) + gaussian(x,0,0.2);
	}
	push_line_plot(&(plot_push_desc){
			.space = &sp,
			.data = data,
			.label = "potential",
			});


	//f64 sample_in[N];
	//for (u32 i = 0; i < N; ++i) {
	//	sample_in[i] = (f64) sp.points[i];
	//}

	//for (u32 i = 0; i < res.component_count; ++i) {
	//	f64 sample_out[N];
	//	settings.basis.sample(res.coeff_count, &res.coeff[i*res.coeff_count], N, sample_out, sample_in);

	//	for (u32 i = 0; i < N; ++i) {
	//		f64 c = fabs(sample_out[i]);
	//		data[i] = c*c;
	//	}
	//	push_line_plot(&(plot_push_desc){
	//			.space = &sp,
	//			.data = data,
	//			.label = plot_snprintf("comp %u", i),
	//			.offset = res.energy[i],
	//			});
	//}
}

void operator(const u32 len, f64 out[static len],
			  f64 in[static len], const u32 component_count,
			  f64 wf[static len*component_count],
			  void* userdata) {
	for (u32 i = 0; i < len; ++i) {
		out[i] = gaussian(in[i], 0.0, 0.2);
		for (u32 j = 0; j < component_count; ++j) {
			f64* g = userdata;
			f64 c = fabs(wf[j*len + i]);
			out[i] += g[j]*c*c;
		}
	}
}

void find_groundstate() {
	if (sp.points == NULL) {
		sp = make_linspace(1, -L/2.0, L/2.0, N);
	}

    struct gp2c_settings settings = {
        .num_basis_functions = ui_num_basis_funcs,
        .max_iterations = ui_max_iterations,
        .error_tol = ui_error_tol,
        //.ho_potential_perturbation = perturbation,
        .gk = gk15,
        .basis = ho_basis,
        .dbgcallback = debug_callback,
        .measure_every = ui_measure_every,
    };

	f64 ga[2] = {-3.0, 1.0};
	f64 gb[2] = { 1.0,-3.0};

    struct gp2c_component components[2] = {
        [0] = {
            .op = operator,
			.userdata = ga,
            //.guess = guess_a
        },
        [1] = {
            .op = operator,
			.userdata = gb,
            //.guess = guess_b
        },
    };

    struct gp2c_result res = gp2c(settings, 2, components);

	/* Plot the result */
	plot_clear();

	for (u32 i = 0; i < N; ++i) {
		f64 x = sp.points[i];
		data[i] = (f32) ho_potential(&x,1,0) + gaussian(x,0,0.2);
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

	for (u32 i = 0; i < res.component_count; ++i) {
		f64 sample_out[N];
		settings.basis.sample(res.coeff_count, &res.coeff[i*res.coeff_count], N, sample_out, sample_in);

		for (u32 i = 0; i < N; ++i) {
			f64 c = fabs(sample_out[i]);
			data[i] = c*c;
		}
		push_line_plot(&(plot_push_desc){
				.space = &sp,
				.data = data,
				.label = plot_snprintf("comp %u", i),
				.offset = res.energy[i],
				});
	}



}

/*
 * UI
 */

void update() {
	igBegin("settings", 0, 0);
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

		struct component_info* p = comp_info_head;
		int iter = 0;
		static char buf[64];
		while (p) {
			snprintf(buf, 64, "component %d", iter);
			if (igButton(buf, (ImVec2){0,0})) {
				p->being_edited = true;
			}
			p = p->next;
			iter++;
		}

		if (igButton("add component", (ImVec2){0,0})) {
			struct component_info** pp = &comp_info_head;
			while (*pp) {
				pp = &(*pp)->next;
			}

			*pp = malloc(sizeof(struct component_info));
			memset(*pp, 0, sizeof(struct component_info));
		}

		if (igButton("find groundstate", (ImVec2){0,0})) {
			find_groundstate();
		}
	igEnd();

	int index = 0;
	for (struct component_info* p = comp_info_head; p; p = p->next) {
		if (p->being_edited) {
			snprintf(buf, 64, "Editing component %d", index);
			igBegin(buf, 0,0);

				//if (igInputInt("basis funcs.", &p->i, 1, 8, 0)) {
				//}

				if (igButton("close", (ImVec2){0,0})) {
					p->being_edited = false;
				}

				if (igButton("delete", (ImVec2){0,0})) {
					int iter = 0;
					struct component_info* node = comp_info_head;
					struct component_info* prev = NULL;
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
}

int main() {
	sbmf_set_log_callback(log_callback);
	sbmf_init();
	plot_init(800, 600, "explore gp");

	plot_set_update_callback(update);
	plot_update_until_closed();

	plot_shutdown();
	sbmf_shutdown();
	return 0;
}
