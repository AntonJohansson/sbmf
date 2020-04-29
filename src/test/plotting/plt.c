#define GLFW_INCLUDE_NONE

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <pthread.h>

#include <GLFW/glfw3.h>
#include <glad/glad.h>

#include "plt.h"
#include "camera.h"
#include "mat4.h"
#include "gl/gl_error_checking.h"
#include "gl/shader.h"
#include "gl/texture.h"

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 600
#define MAX_CURVES_1D 16

#define MAX_GLDRAWDATA_LEN 32

#define MEM_RENDERGROUP_PREALLOC 64*1024

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
} render_entry_type;

typedef struct {
	render_entry_type type;
	uint32_t entry_byte_size;
} render_entry_header;

typedef struct {
	render_entry_header header;

	float* x;
	float* y;
	uint32_t point_count;
} render_entry_line_plot;




typedef struct {
	uint32_t top, block_size;
	uint8_t* memory;
} render_group;

static inline render_group make_render_group(uint32_t size) {
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

static inline void* render_group_push_impl(render_group* group, uint32_t entry_byte_size, render_entry_type type) {
	assert(group->top + entry_byte_size <= group->block_size);

	render_entry_header* header = (render_entry_header*) (group->memory + group->top);
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
	GLenum drawmode;
} gldrawdata;

///////////////////////////////////////////////////////////////////////////////////////

struct PlotState {
	GLFWwindow* window;

	camera cam;

	bool is_lmb_down;
	bool is_mouse_disabled;

	float mouse_x;
	float mouse_y;

	// 1d plotting
	bool active1d[MAX_CURVES_1D];
	unsigned int lengths1d[MAX_CURVES_1D];
	//GLuint vaos1d[MAX_CURVES_1D];
	//GLuint vboxs1d[MAX_CURVES_1D];
	//GLuint vboys1d[MAX_CURVES_1D];
	GLuint program1d;

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

static void mouse_move_callback(GLFWwindow* window, double x, double y){
	PlotState* state = (PlotState*) glfwGetWindowUserPointer(window);

	static double last_normalized_x = 0.0, last_normalized_y = 0.0;

	double normalized_x = 2.0*x/(double)WINDOW_WIDTH - 1.0;
	double normalized_y = 1.0 - 2.0*y/(double)WINDOW_HEIGHT;

	state->mouse_x = normalized_x;
	state->mouse_y = normalized_y;

	double dx = normalized_x - last_normalized_x;
	double dy = normalized_y - last_normalized_y;

	last_normalized_x = normalized_x;
	last_normalized_y = normalized_y;

	if (state->is_lmb_down){
		// If camera is in "arc balL" mode
		state->cam.theta = fclamp(state->cam.theta - dy, 0.0f, M_PI);
		state->cam.phi = fmod(state->cam.phi - dx, 2.0f*M_PI);

		// If camera is in "panning mode"
		state->cam.tar_x+=state->cam.radius*dx;
		state->cam.tar_y+=state->cam.radius*dy;

		//camera_update_arcball(&state->cam);
		camera_update_pan(&state->cam);
	}
}

static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods){
	PlotState* state = (PlotState*) glfwGetWindowUserPointer(window);

	if (button == GLFW_MOUSE_BUTTON_LEFT) {
		switch(action){
			case GLFW_PRESS:{
				state->is_lmb_down = true;
				state->is_mouse_disabled = true;
				glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
				break;
			} 		
			case GLFW_RELEASE:{
				state->is_lmb_down = false; 	
				if (state->is_mouse_disabled){
					state->is_mouse_disabled = false;
					glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
				}
				break;
			}
		}
	}
}

static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset){
	PlotState* state = (PlotState*) glfwGetWindowUserPointer(window);
	state->cam.radius = fmax(state->cam.radius + -0.5f*yoffset, 0.0f);
	//camera_update_arcball(&state->cam);
	camera_update_pan(&state->cam);
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

	GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "plot", 0, 0);
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

static void* plt_update_func(void* data) {
	PlotState* state = (PlotState*) data;

	{
		pthread_mutex_lock(&state->mutex);
		state->window = setup_window();
		setup_gl_error_callback();

		{
			//glGenVertexArrays(MAX_CURVES_1D, state->vaos1d);

			//glGenBuffers(MAX_CURVES_1D, state->vboxs1d);
			//glGenBuffers(MAX_CURVES_1D, state->vboys1d);

			static const char* vert1d = 
				"#version 410\n"
				"layout(location = 0) in float x;\n"
				"layout(location = 1) in float y;\n"
				"uniform mat4 M;\n"
				"uniform mat4 V;\n"
				"uniform mat4 P;\n"
				"void main() {\n"
				"	gl_Position = V*vec4(x, y, 0, 1);\n"
				"}";

			static const char* frag1d =
				"#version 410\n"
				"out vec4 out_color;\n"
				"uniform vec3 color;\n"
				"void main() {\n"
				"	out_color = vec4(color,1);\n"
				"}";

			ShaderProgram prog = program_make();
			program_attach_shader_from_buffer(&prog, vert1d, GL_VERTEX_SHADER);
			program_attach_shader_from_buffer(&prog, frag1d, GL_FRAGMENT_SHADER);
			program_link(&prog);
			program_bind_fragdata_location(&prog, "out_color");
			state->program1d = prog.handle;

		}

		glfwSetWindowUserPointer(state->window, (void*)state);

		glfwSetScrollCallback(state->window, scroll_callback);
		glfwSetCursorPosCallback(state->window, mouse_move_callback);
		glfwSetMouseButtonCallback(state->window, mouse_button_callback);

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
							glGenVertexArrays(1, &gldata->vao);
							glBindVertexArray(gldata->vao);

							glGenBuffers(1, &gldata->vbo);
							glBindBuffer(GL_ARRAY_BUFFER, gldata->vbo);
							// The data in the entry is stored as xxxx....yyyy... where thera are entry->point_count
							// x's and entry->point_count y's.
							glBufferData(GL_ARRAY_BUFFER, sizeof(GL_FLOAT)*2*entry->point_count, entry->x, GL_STATIC_DRAW);
							glUseProgram(state->program1d);
							program_bind_vertex_attribute((VertexAttribute){0,1, GL_FLOAT, GL_FALSE, sizeof(GL_FLOAT), 0});
							program_bind_vertex_attribute((VertexAttribute){1,1, GL_FLOAT, GL_FALSE, sizeof(GL_FLOAT), entry->point_count*sizeof(GL_FLOAT)});

							gldata->vertex_count = entry->point_count;
							gldata->drawmode = GL_LINE_STRIP;
						} break;
						default:
							assert(0);
					};

					ptr += header->entry_byte_size;
				}

				state->group.top = 0;
			}

			while (accumulator >= timestep) {
				accumulator -= timestep;
				plt_update(state);
			}

		}
		pthread_mutex_unlock(&state->mutex);
	}

	{
		//glDeleteBuffers(MAX_CURVES_1D, state->vboxs1d),
		//glDeleteBuffers(MAX_CURVES_1D, state->vboys1d),
		//glDeleteVertexArrays(MAX_CURVES_1D, state->vaos1d);
	}

	//glfwDestroyWindow(state->window);
	glfwTerminate();


	return NULL;
}





























PlotState* plt_init() {
	void* mem = malloc(sizeof(PlotState));
	assert(mem);
	PlotState* state = (PlotState*) mem;
	memset(state, 0, sizeof(PlotState));

	state->cam = make_camera();
	state->cam.radius = 1.0f;
	state->cam.aspect = WINDOW_WIDTH/(float)WINDOW_HEIGHT;
	state->cam.tar_x = 0;
	state->cam.tar_y = 0;
	state->cam.tar_z = 0;
	camera_update_perspective_projection(&state->cam);
	camera_update_orthographic_projection(&state->cam);
	//camera_update_arcball(&state->cam);
	camera_update_pan(&state->cam);


	state->group = make_render_group(MEM_RENDERGROUP_PREALLOC);

	state->gldata_index = 0;

	state->should_join = false;
	pthread_mutex_init(&state->mutex, 0);
	pthread_create(&state->thread, NULL, plt_update_func, state);

	return state;
}

void plt_shutdown(PlotState* state) {
	pthread_mutex_lock(&state->mutex);
	state->should_join = true;
	pthread_mutex_unlock(&state->mutex);

	pthread_join(state->thread, NULL);
	pthread_mutex_destroy(&state->mutex);

	//glDeleteTextures(1, &grid_texture);
	//free_plane_memory();

	//glDeleteProgram(state->program);

	free(state);
}

void plt_update(PlotState* state) {
	glfwPollEvents();

	{
		glClearColor(0,0,0,1);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		//glm::vec4 plane_color = {0.25,0.25,0.25,1};
		//glm::mat4 model(1);
		//model = glm::scale(model, {plane_scale, plane_scale, plane_scale});
		//glm::mat4 mvp = camera.projection * camera.view * model;

		//glUseProgram(screen_program);
		//bind_texture(grid_texture, 1);

		glUseProgram(state->program1d);

		GLuint uniform;
		//uniform = glGetUniformLocation(screen_program, "M");
		//glUniformMatrix4fv(uniform, 1, GL_FALSE, glm::value_ptr(mvp));
		uniform = glGetUniformLocation(state->program1d, "V");
		glUniformMatrix4fv(uniform, 1, GL_TRUE, state->cam.view);
		uniform = glGetUniformLocation(state->program1d, "P");
		glUniformMatrix4fv(uniform, 1, GL_TRUE, state->cam.projection);

		//glEnable(GL_PRIMITIVE_RESTART);
		//glPrimitiveRestartIndex(65535);
		//glDrawElements(GL_TRIANGLE_STRIP, nelements, GL_UNSIGNED_INT, NULL);
		//glDisable(GL_PRIMITIVE_RESTART);

		uniform = glGetUniformLocation(state->program1d, "color");
		//for (unsigned int i = 0; i < MAX_CURVES_1D; ++i) {
		//	//if (!state->active1d[i])
		//	//	continue;

		//	float col[3] = {0};
		//	hsl_to_rgb(col, i*(360.0f/MAX_CURVES_1D),0.85f,0.5f);
		//	glUniform3fv(uniform, 1, col);

		//	//glBindVertexArray(state->vaos1d[i]);
		//	//glDrawArrays(GL_LINE_STRIP, 0, state->lengths1d[i]);
		//}

		for (uint32_t i = 0; i < state->gldata_index; ++i) {
			float col[3] = {0};
			hsl_to_rgb(col, i*(360.0f/state->gldata_index),0.85f,0.5f);
			glUniform3fv(uniform, 1, col);

			gldrawdata* data = &state->gldata[i];
			glBindVertexArray(data->vao);
			glDrawArrays(data->drawmode, 0, data->vertex_count);
		}
	}

	glfwSwapBuffers(state->window);
}

void plt_wait_on_join(PlotState* state) {
	pthread_join(state->thread, NULL);
}




void plt_clear(PlotState* state) {
	pthread_mutex_lock(&state->mutex);
	state->group.top = 0;
	state->gldata_index = 0;
	pthread_mutex_unlock(&state->mutex);
}

void plt_1d(PlotState* state, float* x, float* y, unsigned int len) {
	pthread_mutex_lock(&state->mutex);

	render_entry_line_plot* entry;
	entry = render_group_push_extra_size(&state->group, render_entry_line_plot, 2*len*sizeof(float));
	entry->x = (float*)(entry + 1);
	entry->y = (float*)((uint8_t*)(entry + 1) + len*sizeof(float));
	entry->point_count = len;

	memcpy(entry->x, x, len*sizeof(float));
	memcpy(entry->y, y, len*sizeof(float));

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

void plt_2d(PlotState* state, float* points, unsigned int len);
void plt_3d(PlotState* state, float* points, unsigned int len);
