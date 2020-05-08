#pragma once

#include "common.h"
#include "fdm.h" // for bandmat

extern void eigvp_bandmat(f64* out_eigval, c64* out_eigvec, bandmat bm);
