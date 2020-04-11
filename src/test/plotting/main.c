#define GLFW_INCLUDE_NONE

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "camera.h"
#include "gl/gl_error_checking.h"
#include "gl/shader.h"
#include "gl/texture.h"

#include <bool.h>

#include <stdlib.h>

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 600
#define INITIAL_SUBDIVISIONS 10

#define clamp(val,min,max)\ 
	fmax(fmin(val, max),min)

typedef struct {
	camera cam;

	bool is_lmb_down;
	bool is_mouse_disabled;

	float mouse_x;
	float mouse_y;

	int nelements;
	GLuint vao, vbo, ebo;
	GLuint program;
	GLuint grid_texture;
} plotdata;

static plotdata plt;

static void glfw_error_callback(int code, const char* desc) {
	fprintf(stderr, "GLFW: %s (%d).\n", desc, code);
}

static GLFWwindow* setup_window() {
	if (!glfwInit()) {
		fprintf(stderr, "GLFW: Failed to initialize!\n");
		exit(EXIT_FAILURE);
	}

	glfwSetErrorCallback(glfw_error_callback);

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

	GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "intract", 0, 0);
	if (!window) {
		fprintf(stderr, "GLFW: Failed to create window!\n");
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	glfwMakeContextCurrent(window);

	if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
		fprintf(stderr, "GLAD: Failed to initialize OpenGL context!\n");
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	// Disable VSYNC
	glfwSwapInterval(0);

	//glFrontFace(GL_CW);
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_BACK);
  glEnable(GL_DEPTH_TEST);
  glDepthMask(GL_TRUE);

	printf("Running GLFW %s", glfwGetVersionString());
	printf("Using OpenGL %d.%d", GLVersion.major, GLVersion.minor);

	return window;
}

static void mouse_move_callback(GLFWwindow* window, double x, double y){
	static double last_normalized_x = 0.0, last_normalized_y = 0.0;

	double normalized_x = 2.0*x/(double)WINDOW_WIDTH - 1.0;
	double normalized_y = 1.0 - 2.0*y/(double)WINDOW_HEIGHT;
	mouse_position = {normalized_x, normalized_y};

	double dx = normalized_x - last_normalized_x;
	double dy = normalized_y - last_normalized_y;

	last_normalized_x = normalized_x;
	last_normalized_y = normalized_y;

	if (plt.is_lmb_down){
		plt.cam.theta = clamp(plt.cam.theta - dy, 0.0f, M_PI);
		plt.cam.phi = fmod(plt.cam.phi - dx, 2.0f*M_PI);
		camera_update_arcball(camera);
	}
}
static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods){
	if(button == GLFW_MOUSE_BUTTON_LEFT){
		switch(action){
			case GLFW_PRESS:{
				plt.is_lmb_down = true;
				plt.mouse_disabled = true;
				glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
				break;
			} 		
			case GLFW_RELEASE:{
				plt.is_lmb_down = false; 	
				if (plt.mouse_disabled){
					plt.mouse_disabled = false;
					glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
				}
				break;
			}
		}
	}
}
static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset){
	plt.cam.radius = fmax(plt.cam.radius + -0.5f*yoffset, 0.0f);
	camera_update_arcball(&plt.cam);
}

static float* plane_vbo_memory;
static GLuint* plane_ebo_memory;
static void create_plane_memory(int max_subdivisions){
	void* vbo_mem= malloc(sizeof(float)*(4*(max_subdivisions+1)*(max_subdivisions+1)));
	void* ebo_mem = malloc(sizeof(GLuint)*(2*(max_subdivisions+1)*(max_subdivisions)+(max_subdivisions-1)));

	assert(vbo_mem != NULL && ebo_mem != NULL);
	plane_vbo_memory = (float*)vbo_mem;
	plane_ebo_memory = (GLuint*)ebo_mem;
}
static void free_plane_memory(){
	free(plane_vbo_memory);
	free(plane_ebo_memory);
}

static GLuint subdivided_plane(int subdivisions, GLuint vbo, GLuint ebo){
	//GLfloat vbo_data[2*(subdivisions+1)*(subdivisions+1)];
	//GLuint ebo_data[2*(subdivisions+1)*subdivisions+(subdivisions-1)];
	size_t vbo_size = 4*(subdivisions+1)*(subdivisions+1);
	size_t ebo_size = 2*(subdivisions+1)*subdivisions + (subdivisions-1);

	float x0 = -0.5f, y0 = -0.5f;
	float dx = 1/(float)subdivisions, dy = 1/(float)subdivisions;

	int idx = 0;
	for(int j = 0; j <= subdivisions; ++j){
		for(int i = 0; i <= subdivisions; ++i){
			idx = 4*(i+j*(subdivisions+1));
			// xy
			plane_vbo_memory[idx] 		= x0 + i*dx;
			plane_vbo_memory[idx+1] 	= y0 + j*dy;
			// uv
			plane_vbo_memory[idx+2] = i/(float)subdivisions;
			plane_vbo_memory[idx+3] = 1.0f-j/(float)subdivisions;

			if(j < subdivisions){
				idx = 2*i+j*(2*(subdivisions+1)+1);
				plane_ebo_memory[idx] 	= i+j*(subdivisions+1);
				plane_ebo_memory[idx+1] = i+(j+1)*(subdivisions+1);
			}
		}
		if(j < subdivisions-1){
			idx = 2*(subdivisions+1)+j*(2*(subdivisions+1)+1);
			plane_ebo_memory[idx] = 65535;
		}
	}

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*vbo_size, plane_vbo_memory, GL_DYNAMIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*ebo_size, plane_ebo_memory, GL_DYNAMIC_DRAW);

	return ebo_size;
}

static void setup_gl() {
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	glGenBuffers(1, &vbo);
	glGenBuffers(1, &ebo);

	nelements = subdivided_plane(INITIAL_SUBDIVISIONS, vbo, ebo);

	screen_program = program_begin();
	program_attach_shader_from_file("res/screen.vert", GL_VERTEX_SHADER);
	program_attach_shader_from_file("res/screen.frag", GL_FRAGMENT_SHADER);
	bind_vertex_attribute({0,2,GL_FLOAT,GL_FALSE,4*sizeof(GL_FLOAT),0});
	bind_vertex_attribute({1,2,GL_FLOAT,GL_FALSE,4*sizeof(GL_FLOAT), 2*sizeof(GL_FLOAT)});
	program_link();
	program_end();
}

int main() {
	create_plane_memory(1024);

	GLFWwindow* window = setup_window();
	setup_gl_error_callback();
	setup_gl();

	grid_texture = create_texture(0, INITIAL_SUBDIVISIONS, INITIAL_SUBDIVISIONS, GL_R32F, GL_FLOAT);

  glfwSetCharCallback(window, ImGui_ImplGlfwGL3_CharCallback);
  glfwSetScrollCallback(window, scroll_callback);
  glfwSetCursorPosCallback(window, mouse_move_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);

	plt.cam = make_camera();
	plt.cam.radius = 2.0f;
	plt.cam.aspect = WINDOW_WIDTH/(float)WINDOW_HEIGHT;
	camera_update_projection(&plt.cam);
	camera_update_arcball(&plt.cam);

	double timestep = 1.0/60.0;

	double frame_start = glfwGetTime(), 
				 frame_time = 0.0,
				 frame_end = 0.0,
				 accumulator = 0.0;

	while (!glfwWindowShouldClose(window)) {
		frame_end = glfwGetTime();
		frame_time = frame_end - frame_start;
		frame_start = frame_end;

		accumulator += frame_time;
		while(accumulator >= timestep){
			// update simulation
			accumulator -= timestep;
			glfwPollEvents();

			{
				glClearColor(0,0,0,1);
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

				glm::vec4 plane_color = {0.25,0.25,0.25,1};
				glm::mat4 model(1);
				model = glm::scale(model, {plane_scale, plane_scale, plane_scale});
				glm::mat4 mvp = camera.projection * camera.view * model;

				glUseProgram(screen_program);
				bind_texture(grid_texture, 0);

				GLuint uniform;
				uniform = glGetUniformLocation(screen_program, "mvp");
				glUniformMatrix4fv(uniform, 1, GL_FALSE, glm::value_ptr(mvp));
				uniform = glGetUniformLocation(screen_program, "color");
				glUniform4fv(uniform, 1, glm::value_ptr(plane_color));

				glEnable(GL_PRIMITIVE_RESTART);
				glPrimitiveRestartIndex(65535);
				glDrawElements(GL_TRIANGLE_STRIP, nelements, GL_UNSIGNED_INT, NULL);
				glDisable(GL_PRIMITIVE_RESTART);
			}

			glfwSwapBuffers(window);
		}
	}

	glDeleteTextures(1, &grid_texture);
	free_plane_memory();
	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}
