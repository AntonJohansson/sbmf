#include <sbmf/sbmf.h>
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

void update() {
	igBegin("settings", 0, 0);

	static int num_basis_funcs = 1;

	if (igInputInt("basis funcs.", &num_basis_funcs, 1, 8, 0)) {
		if (num_basis_funcs < 1)
			num_basis_funcs = 1;
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
	igEnd();

	int index = 0;
	for (struct component_info* p = comp_info_head; p; p = p->next) {
		if (p->being_edited) {
			snprintf(buf, 64, "Editing component %d", index);
			igBegin(buf, 0,0);

			if (igInputInt("basis funcs.", &p->i, 1, 8, 0)) {

			}

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
