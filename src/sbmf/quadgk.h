#include "common/common.h"

typedef f64 integrand(f64);

struct segment;
typedef struct segment segment;

struct integrationinfo;
typedef struct integrationinfo integrationinfo;

struct integrationinfo {
	i32 order;

	f64 abs_tol;
	f64 rel_tol;

	i32 max_evals;
	i32 num_evals;

	f64 integral_estimate;
	f64 error_estimate;
};

extern void quadgk(integrand* f, f64 start, f64 end, integrationinfo info);
