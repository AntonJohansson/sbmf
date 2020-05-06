#include "quadgk.h"
#include "prioqueue.h"
#include <string.h> // memcpy
#include <stdio.h> // printf

const f64 xd7[8] = {
	-9.9145537112081263920685469752598e-01,
	-9.4910791234275852452618968404809e-01,
	-8.6486442335976907278971278864098e-01,
	-7.415311855993944398638647732811e-01,
	-5.8608723546769113029414483825842e-01,
	-4.0584515137739716690660641207707e-01,
	-2.0778495500789846760068940377309e-01,
	0.0
};

const f64 wd7[8] = {
	2.2935322010529224963732008059913e-02,
	6.3092092629978553290700663189093e-02,
	1.0479001032225018383987632254189e-01,
	1.4065325971552591874518959051021e-01,
	1.6900472663926790282658342659795e-01,
	1.9035057806478540991325640242055e-01,
	2.0443294007529889241416199923466e-01,
	2.0948214108472782801299917489173e-01
};

const f64 gwd7[4] = {
	1.2948496616886969327061143267787e-01,
	2.797053914892766679014677714229e-01,
	3.8183005050511894495036977548818e-01,
	4.1795918367346938775510204081658e-01
};


struct segment {
	f64 integral_estimate;
	f64 error_estimate;

	f64 start;
	f64 end;
};

f64 norm(f64 r) {
	return r*r;
}

segment evaluate_rule(integrand* f, f64 start, f64 end) {
	static i32 xsize = 8;
	static const f64* x = xd7;
	static const f64* w = wd7;

	static i32 gsize = 4;
	static const f64* gw = gwd7;

	f64 midpoint = 0.5 * (end - start);


	// n1 is:
	// 	0 -- if even order
	// 	1 -- if odd order
	i32 sample_count = xsize;
	i32 n1 = 1 - (sample_count & 1);

	f64 fg, fk;
	f64 Ig = 0, Ik = 0;

	i32 size = gsize - n1;
	for (i32 i = 0; i < size; ++i) {
		fg = f(start + (1+x[2*i+1])*midpoint) + f(start + (1-x[2*i+1])*midpoint);
		fk = f(start + (1+x[2*i-1+1])*midpoint) + f(start + (1-x[2*i-1+1])*midpoint);
		Ig += fg * gw[i];
		Ik += fg * w[2*i+1] + fk * w[2*i-1+1];
	}

	// apparently last point has to be handled differently depending on
	// wether or not the order of integration is odd/even.
	if (n1 == 0) {
		Ik += f(start + midpoint) * w[7]; // 7 == end
	} else {
		f64 f0 = f(start + midpoint);
		Ig += f0 * gw[3]; // 3 == end
		Ik += f0 * w[7] + (f(start + (1+x[7-1])*midpoint) + f(start + (1-x[7-1])*midpoint)) * w[7-1]; // 7-1 == end-1
	}

	f64 Ik_s = Ik*midpoint;
	f64 Ig_s = Ig*midpoint;
	f64 error = norm(Ik_s - Ig_s);

	// add nan/inf check on error
	
	return (segment){Ik_s, error, start, end};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

static bool cmpsegments(void* a, void* b) {
	return (((segment*)a)->error_estimate > ((segment*)b)->error_estimate);
}

void hadapt(integrand* f, f64 start, f64 end, integrationinfo info, segment s) {
	prioqueue* pq = prioqueue_new(32, sizeof(segment), cmpsegments);
	prioqueue_push(pq, &s);

	segment* seg = 0;
	while (	info.error_estimate > info.abs_tol && 
					info.error_estimate > info.rel_tol*norm(info.integral_estimate) && 
					info.num_evals < info.max_evals) {
		seg = (segment*)prioqueue_top(pq);
		prioqueue_pop(pq);

		f64 midpoint = 0.5 * (seg->start + seg->end);
		segment s1 = evaluate_rule(f, seg->start, midpoint);
		segment s2 = evaluate_rule(f, midpoint, seg->end);

		info.integral_estimate = (info.integral_estimate - seg->integral_estimate) + s1.integral_estimate + s2.integral_estimate;
		info.error_estimate = (info.error_estimate - seg->error_estimate) + s1.error_estimate + s2.error_estimate;

		info.num_evals += 4*info.order+2;

		prioqueue_push(pq, &s1);
		prioqueue_push(pq, &s2);

		printf("I: %e E: %e -- (%d/%d)\n", info.integral_estimate, info.error_estimate, info.num_evals, info.max_evals);
	}

	printf("I: %e E: %e\n", info.integral_estimate, info.error_estimate);
	prioqueue_free(pq);
}

void quadgk(integrand* f, f64 start, f64 end, integrationinfo info) {
	segment s = evaluate_rule(f, start, end);
	info.num_evals = 4*info.order + 2;

	info.integral_estimate = s.integral_estimate;
	info.error_estimate = s.error_estimate;

	if (s.error_estimate <= info.abs_tol || 
			s.error_estimate <= info.rel_tol*norm(s.integral_estimate) || 
			info.num_evals >= info.max_evals) {
		printf("I: %e E: %e -- (%d/%d)\n", s.integral_estimate, s.error_estimate, info.num_evals, info.max_evals);
		return;
	}

	hadapt(f, start,end, info, s);
}

/////////////////////////////////////////////////////////////////////////////
