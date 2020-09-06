#pragma once

#include "quadgk.h"

typedef void integrand_vec(f64*,f64*,u32,void*);

integration_result quadgk_vec(integrand_vec* f, f64 start, f64 end, integration_settings settings);
