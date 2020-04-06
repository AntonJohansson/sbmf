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
	int_t dimensions;
	int_t total_pointcount;

	void* memory;
	real_t* mins;
	real_t* maxs;
	int_t* pointcounts;
	real_t* deltas;

	real_t* points;
} grid;

static inline grid generate_grid(int_t dimensions, real_t mins[], real_t maxs[], int_t pointcounts[]) {
	int_t total_pointcount = 1;
	for (int_t i = 0; i < dimensions; ++i)
		total_pointcount *= pointcounts[i];

	int_t size_deltas  = dimensions*sizeof(real_t);
	int_t size_points  = dimensions*total_pointcount*sizeof(real_t);

	void* mem = malloc(dimensions*(2*sizeof(real_t) + sizeof(int_t)) + size_deltas + size_points);
	grid g = {
		.total_pointcount = total_pointcount,
		.dimensions = dimensions,
		.memory = mem,
		.mins 				= (real_t*)(mem),
		.maxs 				= (real_t*)((char*)mem +   dimensions*sizeof(real_t)),
		.pointcounts 	=  (int_t*)((char*)mem + 2*dimensions*sizeof(real_t)),
		.deltas 			= (real_t*)((char*)mem +   dimensions*(2*sizeof(real_t) + sizeof(int_t))),
		.points 			= (real_t*)((char*)mem +   dimensions*(2*sizeof(real_t) + sizeof(int_t)) + size_deltas),
	};

	for (int_t i = 0; i < dimensions; ++i) {
		g.mins[i] = mins[i];
		g.maxs[i] = maxs[i];
		g.pointcounts[i] = pointcounts[i];
		g.deltas[i] = (maxs[i] - mins[i])/pointcounts[i];
	}

	int_t indices[dimensions];
	memset(indices, 0, dimensions*sizeof(int_t));

	for (int_t i = 0; i < total_pointcount; ++i) {
		for (int_t j = 0; j < dimensions; ++j) {
			g.points[dimensions*i + j] = mins[j] + indices[j]*g.deltas[j];
		}

		indices[0] = (indices[0]+1) % pointcounts[0];
		for (int_t j = 0; j < dimensions-1; ++j)
			if (indices[j] % pointcounts[j] == 0)
				indices[j+1] = (indices[j+1]+1) % pointcounts[j+1];
			else
				break;
	}

	return g;
}

static inline grid mimic_grid(grid base) {
	int_t size_deltas  = base.dimensions*sizeof(real_t);
	int_t size_points  = base.dimensions*base.total_pointcount*sizeof(real_t);

	void* mem = malloc(base.dimensions*(2*sizeof(real_t) + sizeof(int_t)) + size_deltas + size_points);
	grid g = {
		.dimensions = base.dimensions,
		.total_pointcount = base.total_pointcount,
		.memory = mem,
		.mins 				= (real_t*)(mem),
		.maxs 				= (real_t*)((char*)mem +   base.dimensions*sizeof(real_t)),
		.pointcounts 	=  (int_t*)((char*)mem + 2*base.dimensions*sizeof(real_t)),
		.deltas 			= (real_t*)((char*)mem +   base.dimensions*(2*sizeof(real_t) + sizeof(int_t))),
		.points 			= (real_t*)((char*)mem +   base.dimensions*(2*sizeof(real_t) + sizeof(int_t)) + size_deltas),
	};

	memcpy(g.mins, base.mins, base.dimensions*(2*sizeof(real_t) + sizeof(int_t)) + size_deltas);

	return g;
}

static inline void free_grid(grid g) {
	free(g.memory);
}

#define FOREACH_GRIDPOINT(grid, i) \
		for (int_t i = 0; i < grid.total_pointcount; i += grid.dimensions)
