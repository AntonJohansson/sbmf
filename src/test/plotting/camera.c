#include "camera.h"
#include <math.h>

void camera_update_arcball(camera* cam) {
	//auto view = glm::mat4(1.0f);
	//view = glm::translate(view, camera.target);
	//view = glm::rotate(view, camera.phi + glm::half_pi<float>(), glm::vec3(0,0,1));
	//view = glm::rotate(view, camera.theta, glm::vec3(1,0,0));
	//view = glm::translate(view, glm::vec3(0,0,camera.radius));

	//camera.position = glm::vec3(view[3]);
	//camera.direction = glm::vec3(0,0,camera.radius);

	//camera.view = glm::inverse(view);

	mat4 view = {0};
	m4identity(&view);
	m4translate(&view, view, cam->tar_x,cam->tar_y,cam->tar_z);
	m4rotate(&view, view, cam->phi + M_PI_2, 0,0,1);
	m4rotate(&view, view, cam->theta, 1,0,0);
	m4translate(&view, view, 0,0,cam->radius);

	// get poisition from translation (last) column
	cam->pos_x = view[3];
	cam->pos_y = view[7];
	cam->pos_z = view[11];

	cam->dir_x = 0;
	cam->dir_y = 0;
	cam->dir_z = cam->radius;

	m4inverse(&cam->view, view);
}

void camera_update_pan(camera* cam) {
	mat4 view = {0};

	m4identity(&view);

	float scale = 1.0/cam->radius;
	m4scale(&view, view, scale, scale, scale);

	m4translate(&view, view, cam->tar_x, cam->tar_y, cam->tar_z);

	m4copy(&cam->view, view);
}

void camera_update_perspective_projection(camera* cam) {
	//camera.projection = glm::perspective(glm::radians(camera.fov), camera.aspect_ratio, camera.near, camera.far);
	mat4* ptr = &cam->projection;
	m4perspective(ptr, cam->fov, cam->aspect, cam->near, cam->far);
}

void camera_update_orthographic_projection(camera* cam) {
	mat4* ptr = &cam->projection;
	m4orthographic(ptr, cam->fov, cam->aspect, cam->near, cam->far);
}

