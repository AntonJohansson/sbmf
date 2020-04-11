#pragma once

#include <math.h>
#include <lapacke.h>

typedef float mat4[16];

#define M4(	m11, m12, m13, m14,						\
						m21, m22, m23, m24,						\
						m31, m32, m33, m34,						\
						m41, m42, m43, m44 ) 					\
	(mat4){	 																\
		m11, m12, m13, m14, 									\
		m21, m22, m23, m24, 									\
		m31, m32, m33, m34, 									\
		m41, m42, m43, m44 										\
	}

static inline void m4copy(mat4* ans, mat4 m) {
	*ans = M4(
			m[0], m[1], m[2], m[3],
			m[4], m[5], m[6], m[7],
			m[8], m[9], m[10], m[11],
			m[12], m[13], m[14], m[15],
			);
}

static inline void m4identity(mat4* ans) {
	*ans = M4(
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1
			);
}

static inline void m4mul(mat4* ans, mat4 a, mat4 b) {
	*ans = M4(
			(a[0] *b[0] + a[1] *b[4]	+ a[2] *b[8] + a[3] *b[12]),  (a[0] *b[1] + a[1] *b[5] + a[2] *b[9] + a[3] *b[13]),  (a[0] *b[2] + a[1] *b[6] + a[2] *b[10] + a[3] *b[14]),  (a[0] *b[3] + a[1] *b[7] + a[2] *b[11] + a[3] *b[15]),
			(a[4] *b[0] + a[5] *b[4]	+ a[6] *b[8] + a[7] *b[12]),  (a[4] *b[1] + a[5] *b[5] + a[6] *b[9] + a[7] *b[13]),  (a[4] *b[2] + a[5] *b[6] + a[6] *b[10] + a[7] *b[14]),  (a[4] *b[3] + a[5] *b[7] + a[6] *b[11] + a[7] *b[15]),
			(a[8] *b[0] + a[9] *b[4]	+ a[10]*b[8] + a[11]*b[12]),  (a[8] *b[1] + a[9] *b[5] + a[10]*b[9] + a[11]*b[13]),  (a[8] *b[2] + a[9] *b[6] + a[10]*b[10] + a[11]*b[14]),  (a[8] *b[3] + a[9] *b[7] + a[10]*b[11] + a[11]*b[15]),
			(a[12]*b[0] + a[13]*b[4]  + a[14]*b[8] + a[15]*b[12]),  (a[12]*b[1] + a[13]*b[5] + a[14]*b[9] + a[15]*b[13]),  (a[12]*b[2] + a[13]*b[6] + a[14]*b[10] + a[15]*b[14]),  (a[12]*b[3] + a[13]*b[7] + a[14]*b[11] + a[15]*b[15])
			);
}

static inline void m4translate(mat4* ans, mat4 m, float x, float y, float z) {
	*ans = M4(
			m[0], 	m[1], 	m[2], 	(x*m[0] + y*m[1] + z*m[2]  + m[3]), 
			m[4], 	m[5], 	m[6], 	(x*m[4] + y*m[5] + z*m[6]  + m[7]), 
			m[8], 	m[9], 	m[10],	(x*m[8] + y*m[9] + z*m[10] + m[11]), 
			m[12], 	m[13], 	m[14], 	m[15], 
			);
}

static inline m4rotate(mat4* ans, mat4 m, float angle, float x, float y, float z) {
	float vec_len = sqrt(x*x + y*y + z*z);
	x /= vec_len;
	y /= vec_len;
	z /= vec_len;

	mat4 rotmat = M4(
				cosf(angle)+x*x*(1-cosf(angle)), 			x*y*(1-cosf(angle))-z*sinf(angle), 			x*z*(1-cosf(angle))+y*sinf(angle), 		0,
				y*x*(1-cosf(angle))+z*sinf(angle) 		cosf(angle)+y*y*(1-cosf(angle)), 				y*z*(1-cosf(angle))-x*sinf(angle) 		0,
				z*x*(1-cosf(angle))-y*sinf(angle), 		z*y*(1-cosf(angle))+x*sinf(angle), 			cosf(angle)+z*z*(1-cosf(angle)), 			0,
				0,																		0,																			0,																		1
			);

	m4mul(ans, rotmat, m);
}

// fov, aspect, near, far
static inline m4perspective(mat4* ans, float fov, float aspect, float near, float far) {
	float uh = 1/tanf(fov/2);
	float uw = uh/aspect;

	*ans = M4(
			uw,	0,  0,										 0,
			0, 	uh, 0, 										 0,
			0, 	0,	-far/(far-near), 			-1,
			0, 	0,  -far*near/(far-near),  0
			);
}

static inline m4inverse(mat4* ans, mat4 m) {
	m4copy(ans, m);

	int pivot[4] = {0};
	LAPAKCE_dgetrf(LAPACKE_ROW_MAJOR, 4,4, *ans, 4, pivot);
	LAPACKE_dgetri(LAPACKE_ROW_MAJOR, 4, *ans, 4, pivot);
}
