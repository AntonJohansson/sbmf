#pragma once

#include <sbmf/methods/gp2c.h>

struct gp2c_result gp2c_gsl(struct gp2c_settings settings, const u32 component_count, struct gp2c_component components[static component_count]);
