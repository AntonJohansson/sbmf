#pragma once

#include <math.h>

typedef struct{
	float radius = 0.0f;
	float theta = 0.0f;
	float phi = 0.0f;

	glm::vec3 position{0,0,0};
	glm::vec3 direction{1,0,0};
	glm::vec3 target{0,0,0};

	float fov = 45.0f;
	float aspect_ratio = 16.0f/9.0f;
	float near = 1.0f;
	float far = 40.0f;

	glm::mat4 projection{1};
	glm::mat4 view{1};
} camera;

static inline camera make_camera() {
	return (camera){
		.radius = 0,
		.theta = 0,
		.phi = 0,

		.pos_x = 0,
		.pos_y = 0,
		.pos_z = 0,

		.dir_x = 0,
		.dir_y = 0,
		.dir_z = 0,

		.tar_x = 0,
		.tar_y = 0,
		.tar_z = 0,

		.fov = M_PI_2,
		.aspect = 16.0f/9.0f,
		.near = 1,
		.far = 10,
	};
}

extern void camera_update_arcball(camera* cam);
extern void camera_update_projection(camera* cam);
