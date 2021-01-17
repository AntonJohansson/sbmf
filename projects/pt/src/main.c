#include <sbmf/sbmf.h>
#include <stdio.h>

#define NA 4
#define NB 0

#define GAA (1.0/3.0)
//#define GAA (-2.0/((f64)NA-1))
#define GAB (+1.0/((f64)NB))
#define GBA (+1.0/((f64)NA))
#define GBB (-2.0/((f64)NB-1))

void log_callback(enum sbmf_log_level log_level, const char* msg) {
	printf("%s\n", msg);
}

int main() {
	sbmf_set_log_callback(log_callback);
	sbmf_init();

	f64 g0[] = {
		GAA, GAB,
		GBA, GBB
	};

	i64 occupations[] = {NA,NB};

	struct nlse_result res = nlse_read_from_binary_file("outbin_50");

	{
		//struct pt_result ptres = rayleigh_schroedinger_pt(res, g0, occupations);
		//struct pt_result ptres = rayleigh_schroedinger_pt_rf(res, 0, g0, occupations);
		//struct pt_result ptres = en_pt_rf(res, 0, g0, occupations);
		struct pt_result ptres = en_pt_rf_new(res, 0, g0, occupations);
		printf("E0:          %.15lf\n", ptres.E0);
		printf("E1:          %.15lf\n", ptres.E1);
		printf("E2:          %.15lf\n", ptres.E2);
		printf("E3:          %.15lf\n", ptres.E3);
		printf("E0+E1:       %.15lf\n", ptres.E0+ptres.E1);
		printf("E0+E1+E2:    %.15lf\n", ptres.E0+ptres.E1+ptres.E2);
		printf("E0+E1+E2+E3: %.15lf\n", ptres.E0+ptres.E1+ptres.E2+ptres.E3);
	}

	sbmf_shutdown();
	return 0;
}
