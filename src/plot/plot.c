#define GLFW_INCLUDE_NONE

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <pthread.h>

#include <GLFW/glfw3.h>
#include <glad/glad.h>

#include "plot.h"
#include "camera.h"
#include "mat4.h"
#include "gl/gl_error_checking.h"
#include "gl/shader.h"
#include "gl/texture.h"

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 600

#define MAX_GLDRAWDATA_LEN 32

#define MEM_RENDERGROUP_PREALLOC 1024*1024

#define fclamp(val,min,max) \
	fmax(fmin(val, max),min)

///////////////////////////////////////////////////////////////////////////////////////
// Why does this abstraction need to exist?
//
//	I'm glad you asked, thank you! I want the drawing of the
//	plots to happen in a separate thread, so we can be in the
//	background, and to do this we can either have all gl* calls
//	in that thread, or switch contexts via glfwMakeContextCurrent(...),
//	however the latter results in a PLATFORM_ERROR, which I don't feel
//	like resolving, mostly cause I don't know how :/. So the former
//	wins by default. Yay!
//
//	Drawing operations are batched into entries, stored in a group
//	and the rendering thread will then convert these entries into
//	the equivalent vbos/vaos, and clear the group.

typedef enum {
	ENTRY_TYPE_render_entry_line_plot,
	ENTRY_TYPE_render_entry_surface_plot,
	ENTRY_TYPE_render_entry_pointcloud_plot,
} render_entry_type;

typedef struct {
	render_entry_type type;
	u32 entry_byte_size;
} render_entry_header;

typedef struct {
	render_entry_header header;

	f32* x;
	f32* y;
	u32 point_count;
} render_entry_line_plot;

typedef struct {
	render_entry_header header;

	f32* x;
	f32* y;
	f32* z;
	u32 point_count;
} render_entry_pointcloud_plot;

typedef struct {
	render_entry_header header;

	f32* x;
	f32* y;
	f32* z;
	u32 point_count;
} render_entry_surface_plot;

typedef struct {
	u32 top, block_size;
	u8* memory;
} render_group;

static inline render_group make_render_group(u32 size) {
	return (render_group) {
		.top = 0,
		.block_size = size,
		.memory = malloc(size),
	};
}

static inline void free_render_group(render_group group) {
	free(group.memory);
}

#define render_group_push(group, type) \
	(type*) render_group_push_impl(group, sizeof(type), ENTRY_TYPE_##type)

#define render_group_push_extra_size(group, type, extra_size) \
	(type*) render_group_push_impl(group, sizeof(type) + extra_size, ENTRY_TYPE_##type)

static inline void* render_group_push_impl(render_group* group,
		uint32_t entry_byte_size, render_entry_type type) {
	assert(group->top + entry_byte_size <= group->block_size);

	render_entry_header* header = (render_entry_header*) (group->memory +
																group->top);
	group->top += entry_byte_size;

	header->type = type;
	header->entry_byte_size = entry_byte_size;

	return header;
}
///////////////////////////////////////////////////////////////////////////////////////

typedef struct {
	uint32_t vertex_count;
	GLuint vao;
	GLuint vbo;

	bool has_ebo;
	GLuint ebo;

	bool has_texture;
	GLuint texture;

	GLuint program;
	GLenum drawmode;
	bool active;
} gldrawdata;

///////////////////////////////////////////////////////////////////////////////////////

struct plotstate {
	GLFWwindow* window;

	camera cam;

	bool is_lmb_down;
	bool is_mouse_disabled;

	float mouse_x;
	float mouse_y;

	float scale_x;
	float scale_y;
	float scale_z;

	render_group group;

	uint32_t gldata_index;
	gldrawdata gldata[MAX_GLDRAWDATA_LEN];

	// Threading data
	pthread_t thread;
	pthread_mutex_t mutex;
	bool should_join;
};

static float hsl_f(float h, float s, float l, int n) {
	float k = fmod((n + h/30), 12);
	float a = s*fmin(l, 1-l);
	return l - a*fmax(fmin(fmin(k-3, 9-k), 1), -1);
}

static void hsl_to_rgb(float outcol[3], float h, float s, float l) {
	outcol[0] = hsl_f(h,s,l, 0);
	outcol[1] = hsl_f(h,s,l, 8);
	outcol[2] = hsl_f(h,s,l, 4);
}

static void glfw_error_callback(int code, const char* desc) {
	fprintf(stderr, "GLFW: %s (%d).\n", desc, code);
}

static void mouse_move_callback(GLFWwindow* window, double x, double y) {
	plotstate* state = (plotstate*) glfwGetWindowUserPointer(window);

	static double last_normalized_x = 0.0, last_normalized_y = 0.0;

	double normalized_x = 2.0*x/(double)WINDOW_WIDTH - 1.0;
	double normalized_y = 1.0 - 2.0*y/(double)WINDOW_HEIGHT;

	state->mouse_x = normalized_x;
	state->mouse_y = normalized_y;

	double dx = normalized_x - last_normalized_x;
	double dy = normalized_y - last_normalized_y;

	last_normalized_x = normalized_x;
	last_normalized_y = normalized_y;

	if (state->is_lmb_down) {
		switch (state->cam.mode) {
			case CAM_PAN: {
				state->cam.tar_x += state->cam.radius*dx;
				state->cam.tar_y += state->cam.radius*dy;
				camera_update_pan(&state->cam);
			}
			break;
			case CAM_ARCBALL: {
				state->cam.theta = fclamp(state->cam.theta - dy, 0.0f, M_PI);
				state->cam.phi = fmod(state->cam.phi - dx, 2.0f*M_PI);
				camera_update_arcball(&state->cam);
			}
			break;
			default:
				assert(0);
		};

	}
}

static void mouse_button_callback(GLFWwindow* window, int button, int action,
																	int mods) {
	plotstate* state = (plotstate*) glfwGetWindowUserPointer(window);

	if (button == GLFW_MOUSE_BUTTON_LEFT) {
		switch(action) {
			case GLFW_PRESS: {
				state->is_lmb_down = true;
				state->is_mouse_disabled = true;
				glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
				break;
			}
			case GLFW_RELEASE: {
				state->is_lmb_down = false;
				if (state->is_mouse_disabled) {
					state->is_mouse_disabled = false;
					glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
				}
				break;
			}
		}
	} else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
		if (state->cam.mode == CAM_PAN)
			camera_set_mode(&state->cam, CAM_ARCBALL);
		else if (state->cam.mode == CAM_ARCBALL)
			camera_set_mode(&state->cam, CAM_PAN);
	}
}

static void scroll_callback(GLFWwindow* window, double xoffset,
														double yoffset) {
	plotstate* state = (plotstate*) glfwGetWindowUserPointer(window);
	state->cam.radius = fmax(state->cam.radius + -0.5f*yoffset, 0.1f);

	switch (state->cam.mode) {
		case CAM_PAN: {
			camera_update_pan(&state->cam);
		}
		break;
		case CAM_ARCBALL: {
			camera_update_arcball(&state->cam);
		}
		break;
		default:
			assert(0);
	};
}

static void  key_callback(GLFWwindow* window, i32 key, i32 scancode, i32 action,
													i32 mods) {
	plotstate* state = (plotstate*) glfwGetWindowUserPointer(window);
	if (action == GLFW_PRESS && key >= GLFW_KEY_0 && key <= GLFW_KEY_9) {
		plot_toggle_active(state, key - GLFW_KEY_0);
	}
	if (key == GLFW_KEY_U) {
		state->scale_x += 0.1f;
	} else if (key == GLFW_KEY_J) {
		state->scale_x -= 0.1f;
	} else if (key == GLFW_KEY_I) {
		state->scale_y += 0.1f;
	} else if (key == GLFW_KEY_K) {
		state->scale_y -= 0.1f;
	} else if (key == GLFW_KEY_O) {
		state->scale_z += 0.1f;
	} else if (key == GLFW_KEY_L) {
		state->scale_z -= 0.1f;
	}
	if (action == GLFW_PRESS && key == GLFW_KEY_R) {
		state->scale_x = 1.0f;
		state->scale_y = 1.0f;
		state->scale_z = 1.0f;
	}
}

//static float* plane_vbo_memory;
//static GLuint* plane_ebo_memory;
//static void create_plane_memory(int max_subdivisions){
//	void* vbo_mem= malloc(sizeof(float)*(4*(max_subdivisions+1)*(max_subdivisions+1)));
//	void* ebo_mem = malloc(sizeof(GLuint)*(2*(max_subdivisions+1)*(max_subdivisions)+(max_subdivisions-1)));
//
//	assert(vbo_mem != NULL && ebo_mem != NULL);
//	plane_vbo_memory = (float*)vbo_mem;
//	plane_ebo_memory = (GLuint*)ebo_mem;
//}
//static void free_plane_memory(){
//	free(plane_vbo_memory);
//	free(plane_ebo_memory);
//}
//
//static GLuint subdivided_plane(int subdivisions, GLuint vbo, GLuint ebo){
//	//GLfloat vbo_data[2*(subdivisions+1)*(subdivisions+1)];
//	//GLuint ebo_data[2*(subdivisions+1)*subdivisions+(subdivisions-1)];
//	size_t vbo_size = 4*(subdivisions+1)*(subdivisions+1);
//	size_t ebo_size = 2*(subdivisions+1)*subdivisions + (subdivisions-1);
//
//	float x0 = -0.5f, y0 = -0.5f;
//	float dx = 1/(float)subdivisions, dy = 1/(float)subdivisions;
//
//	int idx = 0;
//	for(int j = 0; j <= subdivisions; ++j){
//		for(int i = 0; i <= subdivisions; ++i){
//			idx = 4*(i+j*(subdivisions+1));
//			// xy
//			plane_vbo_memory[idx] 		= x0 + i*dx;
//			plane_vbo_memory[idx+1] 	= y0 + j*dy;
//			// uv
//			plane_vbo_memory[idx+2] = i/(float)subdivisions;
//			plane_vbo_memory[idx+3] = 1.0f-j/(float)subdivisions;
//
//			if(j < subdivisions){
//				idx = 2*i+j*(2*(subdivisions+1)+1);
//				plane_ebo_memory[idx] 	= i+j*(subdivisions+1);
//				plane_ebo_memory[idx+1] = i+(j+1)*(subdivisions+1);
//			}
//		}
//		if(j < subdivisions-1){
//			idx = 2*(subdivisions+1)+j*(2*(subdivisions+1)+1);
//			plane_ebo_memory[idx] = 65535;
//		}
//	}
//
//	glBindBuffer(GL_ARRAY_BUFFER, vbo);
//	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*vbo_size, plane_vbo_memory, GL_DYNAMIC_DRAW);
//
//	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
//	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*ebo_size, plane_ebo_memory, GL_DYNAMIC_DRAW);
//
//	return ebo_size;
//}
//

static GLFWwindow* setup_window() {
	glfwSetErrorCallback(glfw_error_callback);

	if (!glfwInit()) {
		fprintf(stderr, "GLFW: Failed to initialize!\n");
		return 0;
	}

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

	//glfwWindowHint(GLFW_SAMPLES, 4);

	GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "plot", 0,
																				0);
	if (!window) {
		fprintf(stderr, "GLFW: Failed to create window!\n");
		glfwTerminate();
		return 0;
	}

	glfwMakeContextCurrent(window);

	if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
		fprintf(stderr, "GLAD: Failed to initialize OpenGL context!\n");
		glfwTerminate();
		return 0;
	}

	// Disable VSYNC
	glfwSwapInterval(0);

	//glFrontFace(GL_CW);
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_BACK);
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);

	glEnable(GL_MULTISAMPLE);

	//printf("Running GLFW %s\n", glfwGetVersionString());
	//printf("Using OpenGL %d.%d\n", GLVersion.major, GLVersion.minor);

	return window;
}

static void* plot_update_func(void* data) {
	plotstate* state = (plotstate*) data;

	static GLuint program2d;
	static GLuint program3d;

	{
		pthread_mutex_lock(&state->mutex);
		state->window = setup_window();
		setup_gl_error_callback();

		{
			//glGenVertexArrays(MAX_CURVES_1D, state->vaos1d);

			//glGenBuffers(MAX_CURVES_1D, state->vboxs1d);
			//glGenBuffers(MAX_CURVES_1D, state->vboys1d);

			static const char* vert2d =
				"#version 410\n"
				"layout(location = 0) in float x;\n"
				"layout(location = 1) in float y;\n"
				"uniform mat4 M;\n"
				"uniform mat4 V;\n"
				"uniform mat4 P;\n"
				"uniform float scale_x;\n"
				"uniform float scale_y;\n"
				"uniform float scale_z;\n"
				"void main() {\n"
				"	gl_Position = P*V*vec4(x, y, 0, 1);\n"
				"}";

			static const char* frag2d =
				"#version 410\n"
				"out vec4 out_color;\n"
				"uniform vec3 color;\n"
				"void main() {\n"
				"	out_color = vec4(color,1);\n"
				"}";

			ShaderProgram prog;
			prog = program_make();
			program_attach_shader_from_buffer(&prog, vert2d, GL_VERTEX_SHADER);
			program_attach_shader_from_buffer(&prog, frag2d, GL_FRAGMENT_SHADER);
			program_link(&prog);
			program_bind_fragdata_location(&prog, "out_color");
			program2d = prog.handle;

			static const char* vert3d =
				"#version 410\n"
				"layout(location = 0) in float x;\n"
				"layout(location = 1) in float y;\n"
				"layout(location = 2) in float z;\n"
				"uniform mat4 M;\n"
				"uniform mat4 V;\n"
				"uniform mat4 P;\n"
				"uniform float scale_x;\n"
				"uniform float scale_y;\n"
				"uniform float scale_z;\n"
				"void main() {\n"
				"	gl_Position = P*V*vec4(scale_x*x, scale_y*y, scale_z*z, 1);\n"
				"}";

			static const char* frag3d =
				"#version 410\n"
				"out vec4 out_color;\n"
				"uniform vec3 color;\n"
				"void main() {\n"
				"	out_color = vec4(color,1);\n"
				"}";

			prog = program_make();
			program_attach_shader_from_buffer(&prog, vert3d, GL_VERTEX_SHADER);
			program_attach_shader_from_buffer(&prog, frag3d, GL_FRAGMENT_SHADER);
			program_link(&prog);
			program_bind_fragdata_location(&prog, "out_color");
			program3d = prog.handle;
		}

		glfwSetWindowUserPointer(state->window, (void*)state);

		glfwSetScrollCallback(state->window, scroll_callback);
		glfwSetCursorPosCallback(state->window, mouse_move_callback);
		glfwSetMouseButtonCallback(state->window, mouse_button_callback);
		glfwSetKeyCallback(state->window, key_callback);

		pthread_mutex_unlock(&state->mutex);
	}

	double timestep = 1.0/60.0;

	double frame_start = glfwGetTime(),
				 frame_time = 0.0,
				 frame_end = 0.0,
				 accumulator = 0.0;

	while (true) {
		frame_end = glfwGetTime();
		frame_time = frame_end - frame_start;
		frame_start = frame_end;

		accumulator += frame_time;

		pthread_mutex_lock(&state->mutex);
		{
			if (glfwWindowShouldClose(state->window) || state->should_join) {
				pthread_mutex_unlock(&state->mutex);
				break;
			}

			// Push queued render entries to GL-representation
			if (state->group.top > 0) {

				uint8_t* ptr = state->group.memory;
				uint32_t iteration = 0;
				while (ptr < state->group.memory + state->group.top && iteration < 100) {
					iteration++;

					render_entry_header* header = (render_entry_header*) ptr;
					switch (header->type) {
						case ENTRY_TYPE_render_entry_line_plot: {
							//glBindVertexArray(state->vaos1d[id]);

							//glBindBuffer(GL_ARRAY_BUFFER, state->vboxs1d[id]);
							//glBufferData(GL_ARRAY_BUFFER, sizeof(float)*len, x, GL_STATIC_DRAW);
							//program_bind_vertex_attribute((VertexAttribute){0,1,GL_FLOAT,GL_FALSE,sizeof(GL_FLOAT),0});
							//glBindBuffer(GL_ARRAY_BUFFER, state->vboys1d[id]);
							//glBufferData(GL_ARRAY_BUFFER, sizeof(float)*len, y, GL_STATIC_DRAW);
							//program_bind_vertex_attribute((VertexAttribute){1,1,GL_FLOAT,GL_FALSE,sizeof(GL_FLOAT),0});

							//state->lengths1d[id] = len;
							//state->active1d[id] = true;

							render_entry_line_plot* entry = (render_entry_line_plot*) header;

							assert(state->gldata_index + 1 < MAX_GLDRAWDATA_LEN);
							gldrawdata* gldata = &state->gldata[state->gldata_index++];
							gldata->active = true;
							glGenVertexArrays(1, &gldata->vao);
							glBindVertexArray(gldata->vao);

							glGenBuffers(1, &gldata->vbo);
							glBindBuffer(GL_ARRAY_BUFFER, gldata->vbo);
							// The data in the entry is stored as xxxx....yyyy... where thera are entry->point_count
							// x's and entry->point_count y's.
							glBufferData(GL_ARRAY_BUFFER, sizeof(GL_FLOAT)*2*entry->point_count, entry->x,
													 GL_STATIC_DRAW);
							glUseProgram(program2d);
							program_bind_vertex_attribute((VertexAttribute) {
								0,1, GL_FLOAT, GL_FALSE, sizeof(GL_FLOAT), 0
							});
							program_bind_vertex_attribute((VertexAttribute) {
								1,1, GL_FLOAT, GL_FALSE, sizeof(GL_FLOAT), entry->point_count* sizeof(GL_FLOAT)
							});

							gldata->vertex_count = entry->point_count;
							gldata->drawmode = GL_LINE_STRIP;
							gldata->program = program2d;
						}
						break;
						case ENTRY_TYPE_render_entry_pointcloud_plot: {
							render_entry_pointcloud_plot* entry = (render_entry_pointcloud_plot*) header;

							assert(state->gldata_index + 1 < MAX_GLDRAWDATA_LEN);
							gldrawdata* gldata = &state->gldata[state->gldata_index++];
							gldata->active = true;

							glGenVertexArrays(1, &gldata->vao);
							glBindVertexArray(gldata->vao);

							glGenBuffers(1, &gldata->vbo);
							glBindBuffer(GL_ARRAY_BUFFER, gldata->vbo);
							// The data in the entry is stored as xxxx....yyyy...zzzzzz.... where thera are entry->point_count
							// x's and entry->point_count y's.
							glBufferData(GL_ARRAY_BUFFER, sizeof(GL_FLOAT)*3*entry->point_count, entry->x,
													 GL_STATIC_DRAW);
							glUseProgram(program3d);
							program_bind_vertex_attribute((VertexAttribute) {
								0,1, GL_FLOAT, GL_FALSE, sizeof(GL_FLOAT), 0
							});
							program_bind_vertex_attribute((VertexAttribute) {
								1,1, GL_FLOAT, GL_FALSE, sizeof(GL_FLOAT), entry->point_count* sizeof(GL_FLOAT)
							});
							program_bind_vertex_attribute((VertexAttribute) {
								2,1, GL_FLOAT, GL_FALSE, sizeof(GL_FLOAT),
								2*entry->point_count* sizeof(GL_FLOAT)
							});

							gldata->vertex_count = entry->point_count;
							gldata->drawmode = GL_POINTS;
							gldata->program = program3d;
						}
						break;
						default:
							assert(0);
					};

					ptr += header->entry_byte_size;
				}

				state->group.top = 0;
			}

			while (accumulator >= timestep) {
				accumulator -= timestep;
				plot_update(state);
			}

		}
		pthread_mutex_unlock(&state->mutex);
	}

	{
		//glDeleteBuffers(MAX_CURVES_1D, state->vboxs1d),
		//glDeleteBuffers(MAX_CURVES_1D, state->vboys1d),
		//glDeleteVertexArrays(MAX_CURVES_1D, state->vaos1d);
	}

	glfwTerminate();

	return NULL;
}





























plotstate* plot_init() {
	void* mem = malloc(sizeof(plotstate));
	assert(mem);
	plotstate* state = (plotstate*) mem;
	memset(state, 0, sizeof(plotstate));

	state->cam = make_camera();
	state->cam.radius = 1.0f;
	state->cam.aspect = WINDOW_WIDTH/(float)WINDOW_HEIGHT;
	state->cam.tar_x = 0;
	state->cam.tar_y = 0;
	state->cam.tar_z = 0;
	state->scale_x = 1;
	state->scale_y = 1;
	state->scale_z = 1;

	camera_set_mode(&state->cam, CAM_PAN);

	state->group = make_render_group(MEM_RENDERGROUP_PREALLOC);

	state->gldata_index = 0;

	state->should_join = false;
	pthread_mutex_init(&state->mutex, 0);
	pthread_create(&state->thread, NULL, plot_update_func, state);

	return state;
}

void plot_shutdown(plotstate* state) {
	pthread_mutex_lock(&state->mutex);
	state->should_join = true;
	pthread_mutex_unlock(&state->mutex);

	pthread_join(state->thread, NULL);
	pthread_mutex_destroy(&state->mutex);

	//glDeleteTextures(1, &grid_texture);
	//free_plane_memory();

	//glDeleteProgram(state->program);

	free_render_group(state->group);
	free(state);
}

void plot_update(plotstate* state) {
	glfwPollEvents();

	{
		glClearColor(0,0,0,1);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		for (uint32_t i = 0; i < state->gldata_index; ++i) {
			gldrawdata* data = &state->gldata[i];

			if (!data->active)
				continue;

			glUseProgram(data->program);

			GLuint uniform;
			//uniform = glGetUniformLocation(screen_program, "M");
			//glUniformMatrix4fv(uniform, 1, GL_FALSE, glm::value_ptr(mvp));
			uniform = glGetUniformLocation(data->program, "V");
			glUniformMatrix4fv(uniform, 1, GL_TRUE, state->cam.view);
			uniform = glGetUniformLocation(data->program, "P");
			glUniformMatrix4fv(uniform, 1, GL_TRUE, state->cam.projection);

			uniform = glGetUniformLocation(data->program, "scale_x");
			glUniform1fv(uniform, 1, &state->scale_x);
			uniform = glGetUniformLocation(data->program, "scale_y");
			glUniform1fv(uniform, 1, &state->scale_y);
			uniform = glGetUniformLocation(data->program, "scale_z");
			glUniform1fv(uniform, 1, &state->scale_z);

			uniform = glGetUniformLocation(data->program, "color");
			float col[3] = {0};
			hsl_to_rgb(col, i*(360.0f/state->gldata_index),0.85f,0.5f);
			glUniform3fv(uniform, 1, col);

			glBindVertexArray(data->vao);
			glDrawArrays(data->drawmode, 0, data->vertex_count);
		}
	}

	glfwSwapBuffers(state->window);
}

void plot_wait_on_join(plotstate* state) {
	pthread_join(state->thread, NULL);
}

void plot_toggle_active(plotstate* state, i32 idx) {
	if (idx >= 0 && idx < state->gldata_index) {
		state->gldata[idx].active = !state->gldata[idx].active;
	}
}



void plot_clear(plotstate* state) {
	pthread_mutex_lock(&state->mutex);
	state->group.top = 0;
	state->gldata_index = 0;
	pthread_mutex_unlock(&state->mutex);
}

void plot_1d(plotstate* state, f32* x, f32* y, u32 len) {
	pthread_mutex_lock(&state->mutex);

	render_entry_line_plot* entry;
	entry = render_group_push_extra_size(&state->group, render_entry_line_plot,
																			 2*len*sizeof(f32));
	entry->x = (f32*)(entry + 1);
	entry->y = (f32*)((u8*)(entry + 1) + len*sizeof(f32));
	entry->point_count = len;

	memcpy(entry->x, x, len*sizeof(f32));
	memcpy(entry->y, y, len*sizeof(f32));

	//glfwMakeContextCurrent(state->window);

	//glBindVertexArray(state->vaos1d[id]);

	//glBindBuffer(GL_ARRAY_BUFFER, state->vboxs1d[id]);
	//glBufferData(GL_ARRAY_BUFFER, sizeof(float)*len, x, GL_STATIC_DRAW);
	//program_bind_vertex_attribute((VertexAttribute){0,1,GL_FLOAT,GL_FALSE,sizeof(GL_FLOAT),0});
	//glBindBuffer(GL_ARRAY_BUFFER, state->vboys1d[id]);
	//glBufferData(GL_ARRAY_BUFFER, sizeof(float)*len, y, GL_STATIC_DRAW);
	//program_bind_vertex_attribute((VertexAttribute){1,1,GL_FLOAT,GL_FALSE,sizeof(GL_FLOAT),0});

	//state->lengths1d[id] = len;
	//state->active1d[id] = true;

	pthread_mutex_unlock(&state->mutex);
}

void plot_2d(plotstate* state, f32* x, f32* y, f32* z, u32 len) {
	pthread_mutex_lock(&state->mutex);

	render_entry_pointcloud_plot* entry;
	entry = render_group_push_extra_size(&state->group, render_entry_pointcloud_plot, 3*len*sizeof(f32));
	entry->x = (f32*)(entry + 1);
	entry->y = (f32*)((u8*)(entry + 1) +   len*sizeof(f32));
	entry->z = (f32*)((u8*)(entry + 1) + 2*len*sizeof(f32));
	entry->point_count = len;

	memcpy(entry->x, x, len*sizeof(f32));
	memcpy(entry->y, y, len*sizeof(f32));
	memcpy(entry->z, z, len*sizeof(f32));

	pthread_mutex_unlock(&state->mutex);
}

void plot_3d(plotstate* state, f64* points, u32 len) {}






plotsurface* plot_make_surface(plotstate* state, u32 size) {
	u32 vbo_size = sizeof(f32)*(4*(size)*(size));
	u32 ebo_size = sizeof(u32)*(2*(size)* (size - 1)+(size-2));

	plotsurface* surf = xmalloc(sizeof(plotsurface) + vbo_size + ebo_size);
	surf->size = size;
	surf->vbo_mem = (f32*)(surf+1);
	surf->ebo_mem = (f32')((u8*)(surf->vbo_mem)+vbo_size);


	f32 x0 = -0.5f, y0 = -0.5f;
	f32 dx = 1/(f32)subdivisions, dy = 1/(f32)subdivisions;

	i32 idx = 0;
	for(i32 j = 0; j < size; ++j){
		for(i32 i = 0; i < size; ++i){
			idx = 4*(i+j*(size));
			// xy
			plane_vbo_memory[idx] 		= x0 + i*dx;
			plane_vbo_memory[idx+1] 	= y0 + j*dy;
			// uv
			plane_vbo_memory[idx+2] = i/(f32)(size-1);
			plane_vbo_memory[idx+3] = 1.0f-j/(f32)(size-1);

			if(j < size-1){
				idx = 2*i+j*(2*(size)+1);
				plane_ebo_memory[idx] 	= i+j*(size);
				plane_ebo_memory[idx+1] = i+(j+1)*(size);
			}
		}

		if(j < size-2){
			idx = 2*(size)+j*(2*(size)+1);
			plane_ebo_memory[idx] = 65535;
		}
	}

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*vbo_size, plane_vbo_memory, GL_DYNAMIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*ebo_size, plane_ebo_memory, GL_DYNAMIC_DRAW);

	return surf;
}

void plot_free_surface(plotstate* state, plotsurface* surf) {
	free(surf);
}











void plot_surface(plotstate* state, plotsurface* surf, f32* z) {
	pthread_mutex_lock(&state->mutex);


	render_entry_surface_plot* entry;
	entry = render_group_push_extra_size(&state->group, render_entry_surface_plot, vbo_size + ebo_size);
	pthread_mutex_unlock(&state->mutex);
}
