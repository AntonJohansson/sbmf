#pragma once

#include "common.h"
#include <stdlib.h>
#include <string.h>

// For example, consider the block of lengths 1x2x3 sampled with 128,64,32 points in each dimension, this
// corresponds to the grid with
//
//  	dimension = 3,
//  	lengths = {1, 2, 3},
//  	points_per_length = {128,64,32},
//  	deltas = {1/128, 2/64, 3/32},
//  	points = {...} (128 x 64 x 32 array)

typedef struct {
	real_t min, max;
	int_t pointcount;
} dimension_info;

typedef struct {
	void* memory;

	int_t dimensions;
	dimension_info* infos;

	real_t* deltas;

	int_t total_pointcount;
	real_t* points;
} grid;

static inline grid generate_grid(int_t dimensions, dimension_info infos[]) {
	int_t total_pointcount = 1;
	for (int_t i = 0; i < dimensions; ++i)
		total_pointcount *= infos[i].pointcount;

	int_t size_diminfo = dimensions*sizeof(dimension_info);
	int_t size_deltas  = dimensions*sizeof(real_t);
	int_t size_points  = dimensions*total_pointcount*sizeof(real_t);

	void* mem = malloc(size_diminfo + size_deltas + size_points);
	grid g = {
		.memory = mem,
		.dimensions = dimensions,
		.infos = (dimension_info*)mem,
		.deltas = (real_t*)((char*)mem + size_diminfo),
		.total_pointcount = total_pointcount,
		.points = (real_t*)((char*)mem + size_diminfo + size_deltas),
	};

	for (int_t i = 0; i < dimensions; ++i) {
		g.infos[i] = infos[i];
		g.deltas[i] = (g.infos[i].max - g.infos[i].min)/g.infos[i].pointcount;
	}

	int_t indices[dimensions];
	memset(indices, 0, dimensions*sizeof(int_t));

	for (int_t i = 0; i < total_pointcount; ++i) {
		for (int_t j = 0; j < dimensions; ++j) {
			g.points[dimensions*i + j] = g.infos[j].min + indices[j]*g.deltas[j];
		}

		indices[0] = (indices[0]+1) % g.infos[0].pointcount;
		for (int_t j = 0; j < dimensions-1; ++j)
			if (indices[j] % g.infos[j].pointcount == 0)
				indices[j+1] = (indices[j+1]+1) % g.infos[j+1].pointcount;
			else
				break;
	}

	return g;
}

static inline void free_grid(grid g) {
	free(g.memory);
}
