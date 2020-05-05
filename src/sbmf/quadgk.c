#include "quadgk.h"

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
	static const f64 xk[] = xd7;
	static const f64 w[] = wd7;
	static const f64 gw[] = gwd7;

	f64 midpoint = 0.5 * (end - start);


	// n1 is:
	// 	0 -- if even order
	// 	1 -- if odd order
	i32 sample_count = sizeof(x)/sizeof(x[0]);
	i32 n1 = 1 - (sample_count & 1);

	f64 fg, fk;
	f64 Ig = 0, Ik = 0;

	i32 size = sizeof(gw)/sizeof(gw[0]) - n1;
	for (i32 i = 0; i < size; ++i) {
		fg = f(start + (1+x[2*i+1])*midpoint) + f(start + (1-x[2*i+1])*midpoint);
		fk = f(start + (1+x[2*i-1+1])*midpoint) + f(start + (1-x[2*i-1+1])*midpoint);
		Ig += fg * gw[i];
		Ik += fg * w[2*i+1] + fk * w[2*i-1+1];
	}

	f64 IK_s = Ik*midpoint;
	f64 Ig_s = Ig*midpoint;
	f64 error = norm(Ik_s - Ig_s);

	return (segment){Ik_s, error, start, end};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

struct integrationinfo {
	i32 order;

	f64 abs_tol;
	f64 rel_tol;

	i32 max_evals;
	i32 num_evals;

	f64 integral_estimate;
	f64 error_estimate;
};

void hadapt(integrand* f, f64 start, f64 end, integrationinfo info, segment s) {
}

void quadgk(integrand* f, f64 start, f64 end, integrationinfo info) {
	segment s = evaluate_rule(f, start, end);
	info.numevals = 4*info.order + 2;

	info.integral_estimate = s.integral_estimate;
	info.error_estimate = s.error_estimate;

	if (s.error_estimate <= info.abs_tol || 
			s.error_estimate <= info.rel_tol*norm(s.integral_estimate) || 
			info.numevals >= info.maxevals) {
		printf("I: %e E: -- (%d/%d)\n", s.integral_estimate, s.error_estimate, info.numevals, info.maxevals);
		return;
	}

	hadapt(f, a,b, info, s);
}

/////////////////////////////////////////////////////////////////////////////
