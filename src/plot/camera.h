#pragma once

#include "mat4.h"
#include "common.h"

typedef struct {
	f32 radius;
	f32 theta;
	f32 phi;

	f32 pos_x, pos_y, pos_z;
	f32 dir_x, dir_y, dir_z;
	f32 tar_x, tar_y, tar_z;
	
	f32 fov;
	f32 aspect;
	f32 near, far;

	mat4 projection;
	mat4 view;

	enum {
		CAM_PAN = 0,
		CAM_ARCBALL = 1,
	} mode;
} camera;

static inline camera make_camera() {
	return (camera){
		.radius = 1,
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
		.aspect = 4.0/3.0,
		.near = 0.01,
		.far = 50,

		.mode = CAM_PAN,
	};
}

extern void camera_update_arcball(camera* cam);
extern void camera_update_pan(camera* cam);
extern void camera_update_perspective_projection(camera* cam);
extern void camera_update_orthographic_projection(camera* cam);
extern void camera_set_mode(camera* cam,i32 mode);
