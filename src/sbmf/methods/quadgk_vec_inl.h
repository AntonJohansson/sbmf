#pragma once

#include "quadgk_vec.h"

integration_result quadgk_vec_inl(integrand_vec* f, f64 start, f64 end, integration_settings settings);
