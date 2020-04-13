#pragma once

#include <math.h>
#include "mat4.h"

typedef struct {
	float radius;
	float theta;
	float phi;

	float pos_x, pos_y, pos_z;
	float dir_x, dir_y, dir_z;
	float tar_x, tar_y, tar_z;
	
	float fov;
	float aspect;
	float near, far;

	mat4 projection;
	mat4 view;
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

		.fov = M_PI_4,
		.aspect = 16.0f/9.0f,
		.near = 1,
		.far = 10,
	};
}

extern void camera_update_arcball(camera* cam);
extern void camera_update_pan(camera* cam);
extern void camera_update_perspective_projection(camera* cam);
extern void camera_update_orthographic_projection(camera* cam);
