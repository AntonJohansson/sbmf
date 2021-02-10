static inline void map_to_triangular_index(u32 k, u32 N, u32* m, u32* n) {
	*m = k / N;
	*n = k % N;
	if (*m > *n) {
		*m = N - *m - 0;
		*n = N - *n - 1;
	}
}
/*
 * Holds all information needed to do the PT,
 * easy to pass around. Basicly params to
 * pt_rayleigh_schroedinger
 */
struct pt_settings {
	struct nlse_result* res;
	struct nlse_settings* settings;
	f64* g0;
	i64* particle_count;

	struct eigen_result_real* states;
	const u32 N; /* states to include */
	const u32 L; /* coeff count */
	nlse_operator_func* pert;
	f64* hermite_cache;
};

__host__ __device__
static inline f64 G0(struct pt_settings* pt, u32 A, u32 B) {
	return pt->g0[A*pt->res->component_count + B];
}

__host__ __device__
static inline f64 E(struct pt_settings* pt, u32 A, u32 i) {
	return pt->states[A].eigenvalues[i];
}

__host__ __device__
static inline f64* PHI(struct pt_settings* pt, u32 A, u32 i) {
	return &pt->states[A].eigenvectors[i * pt->L];
}

__host__ __device__
static inline f64 V_closed(struct pt_settings* pt, f64* cache, u32 A, u32 B, u32 i, u32 j, u32 k, u32 l) {
	f64* phi_a = PHI(pt, A, i);
	f64* phi_b = PHI(pt, B, j);
	f64* phi_c = PHI(pt, A, k);
	f64* phi_d = PHI(pt, B, l);

	f64 sum = 0.0;
	for (u32 a = 0; a < pt->L; ++a) {
		for (u32 b = 0; b < pt->L; ++b) {
			for (u32 c = 0; c < pt->L; ++c) {
				for (u32 d = 0; d < pt->L; ++d) {
					f64 L = phi_a[a]*phi_b[b]*phi_c[c]*phi_d[d];//*ho_K(a)*ho_K(b)*ho_K(c)*ho_K(d);
					if (fabs(L) < 1e-10)
						continue;
					f64 integral = cache[index4_cuda(a,b,c,d)];
					sum += L*integral;
				}
			}
		}
	}

	return sum;
}

/*
 * Helper functions since these will be calculated a lot
 */

__host__ __device__
static inline f64 rs_2nd_order_me(struct pt_settings* pt, u32 A, u32 B, u32 i, u32 j) {
	f64 me = 0.0;
	if (A == B) {
		f64 factor = (i == j) ? 1.0/sqrt(2.0) : 1.0;
		me = factor * G0(pt,A,A) * sqrt((f64)(pt->particle_count[A] * (pt->particle_count[A] - 1)))
			* V_closed(pt, pt->hermite_cache, A,A, i,j,0,0);
	} else {
		me = G0(pt,A,B) * sqrt((f64) (pt->particle_count[A] * pt->particle_count[B]))
			* V_closed(pt, pt->hermite_cache, A,B, i,j,0,0);
	}
	return me;
}

__host__ __device__
static inline f64 rs_2nd_order_ediff(struct pt_settings* pt, u32 A, u32 B, u32 i, u32 j) {
	return E(pt,A,0) + E(pt,B,0) - E(pt,A,i) - E(pt,B,j);
}

/*
 * Main function for Rayleigh-Schrodinger perturbation theory
 */

struct pt_result rspt_1comp_cuda(struct nlse_settings settings, struct nlse_result res, u32 component, f64* g0, i64* particle_count) {
	/* order of hamiltonians, that is include all states */
	const u32 states_to_include = res.coeff_count;
	const i64* N = particle_count;
	sbmf_log_info("running 1comp RSPT cuda:\n    components: %u\n    states: %u\n", res.component_count, states_to_include);

	struct eigen_result_real states;
	states = find_eigenpairs_full_real(res.hamiltonian[component]);
	for (u32 j = 0; j < states_to_include; ++j) {
		f64_normalize(&states.eigenvectors[j*res.coeff_count], &states.eigenvectors[j*res.coeff_count], res.coeff_count);
	}

	const u64 hermite_integral_count = size4_cuda(states_to_include-1);
	const u64 hermite_cache_size = sizeof(f64)*hermite_integral_count;
	f64* hermite_cache = (f64*)sbmf_stack_push(hermite_cache_size);
	{
		sbmf_log_info("Precomputing %ld hermite integrals", hermite_integral_count);
		for (u32 i = 0; i < states_to_include; ++i) {
			for (u32 j = i; j < states_to_include; ++j) {
				for (u32 k = j; k < states_to_include; ++k) {
					for (u32 l = k; l < states_to_include; ++l) {
						hermite_cache[index4_cuda(i,j,k,l)] = hermite_integral_4_cuda(i,j,k,l);
					}
				}
			}
		}
	}

	f64* hermite_cache_device;
	cudaMalloc(&hermite_cache_device, hermite_cache_size);
	cudaMemcpy(hermite_cache_device, hermite_cache, hermite_cache_size, cudaMemcpyHostToDevice);

	struct pt_settings pt = {
		.res = &res,
		.settings = &settings,
		.g0 = g0,
		.particle_count = particle_count,
		.states = &states,
		.N = states_to_include,
		.L = res.coeff_count,
		.pert = settings.spatial_pot_perturbation,
		.hermite_cache = hermite_cache,
	};

	/* Zeroth order PT */
	sbmf_log_info("Starting zeroth order PT");
	f64 E0 = 0.0;
	{
		E0 += N[component] * E(&pt, component, 0);
	}
	sbmf_log_info("\tE0: %e", E0);

	/* First order PT */
	sbmf_log_info("Starting first order PT");
	f64 E1 = 0.0;
	{
		/* Handles interaction within component */
		E1 += -0.5 * G0(&pt,component,component) * N[component] * (N[component]-1) * V_closed(&pt, hermite_cache, component,component, 0,0,0,0);
	}
	sbmf_log_info("\tE1: %e", E1);

	const u32 pt2_cache_size = ((states_to_include-1)*(states_to_include))/2;
	f64 pt2_cache[pt2_cache_size];

	/* Assumes i in [0,states_to_include), j in [0,states_to_include) */
#define PT2_CACHE_INDEX(i, j) \
	((i)*states_to_include - (((i)*(i+1))/2) + j)


	/* Second order PT */
	sbmf_log_info("Starting second order PT");
	f64 E2 = 0.0;
	{
		/*
		 * Double substitutions (both excitations within same component),
		 * loop over unique pairs (j,k).
		 */
#pragma omp parallel for reduction(+: E2)
		for (u32 i = 1; i < states_to_include; ++i) {
			for (u32 j = i; j < states_to_include; ++j) {
				f64 me = rs_2nd_order_me(&pt, component,component, i,j);
				f64 Ediff = rs_2nd_order_ediff(&pt, component,component, i,j);

				pt2_cache[PT2_CACHE_INDEX(i-1, j-i)] = me/Ediff;

				E2 += me*me/(Ediff);
			}
		}
	}
	sbmf_log_info("\tE2: %e", E2);

	/* Third order PT */
	sbmf_log_info("Starting third order PT");
	f64 E3 = 0.0;
	{
		f64 E_00_00 = 0;
		{
#pragma omp parallel for reduction(+: E_00_00)
			for (u32 m = 1; m < states_to_include; ++m) {
				for (u32 n = m; n < states_to_include; ++n) {
					f64 tmn = pt2_cache[PT2_CACHE_INDEX(m-1,n-m)];
					E_00_00 += tmn*tmn;
				}
			}

			f64 v_00_00 = G0(&pt,component,component)*V_closed(&pt, hermite_cache, component,component, 0,0,0,0);
			E_00_00 *= v_00_00;
		}
		sbmf_log_info("\t\t00,00: %.10e", E_00_00);

		f64 E_m0_n0 = 0;
		{
			const f64 c_root_2_minus_1 = sqrt(2.0) - 1.0;
			const f64 c_3_minus_2_root_2 = 3.0 - 2.0*sqrt(2.0);

#pragma omp parallel for reduction(+: E_m0_n0)
			for (u32 k = 0; k < ((states_to_include-1)*states_to_include)/2; ++k) {
				u32 m = k/(states_to_include-1);
				u32 n = k%(states_to_include-1);
				if (m > n) {
					m = (states_to_include-1) - m - 0;
					n = (states_to_include-1) - n - 1;
				}
				m += 1;
				n += 1;

				const f64 v_m0_n0 = G0(&pt,component,component)*V_closed(&pt, hermite_cache, component,component, m,0,n,0);

				f64 sum = 0.0;
				for (u32 p = 1; p < states_to_include; ++p) {

					f64 tmp = 0;
					if (p >= m)
						tmp = pt2_cache[PT2_CACHE_INDEX(m-1, p-m)];
					else
						tmp = pt2_cache[PT2_CACHE_INDEX(p-1, m-p)];

					f64 tnp = 0;
					if (p >= n)
						tnp = pt2_cache[PT2_CACHE_INDEX(n-1, p-n)];
					else
						tnp = pt2_cache[PT2_CACHE_INDEX(p-1, n-p)];

					const f64 delta_mp = (m == p) ? 1.0 : 0.0;
					const f64 delta_np = (n == p) ? 1.0 : 0.0;

					const f64 coeff = 1 + c_root_2_minus_1*(delta_mp + delta_np) + c_3_minus_2_root_2*(delta_mp*delta_np);
					sum += coeff * tmp * tnp;
				}

				f64 factor = (m == n) ? 1.0 : 2.0;

				E_m0_n0 += factor * v_m0_n0 * sum;
			}

			E_m0_n0 *= (N[component] - 3);
		}
		sbmf_log_info("\t\tm0,n0: %.10e", E_m0_n0);

		f64 E_mn_pq = 0;
		{
			const f64 c_root_2_minus_2 = sqrt(2.0) - 2.0;
			const f64 c_3_minus_2_root_2 = 3.0 - 2.0*sqrt(2.0);

			const u32 INDS_N = ((states_to_include-1)*states_to_include)/2;

#pragma omp parallel for reduction(+: E_mn_pq)
			for (u32 k = 0; k < INDS_N*(INDS_N+1)/2; ++k) {
				u32 k0, k1;
				map_to_triangular_index(k, INDS_N, &k0, &k1);

				u32 m, n;
				map_to_triangular_index(k0, states_to_include-1, &m, &n);
				m += 1; n += 1;

				u32 p, q;
				map_to_triangular_index(k1, states_to_include-1, &p, &q);
				p += 1; q += 1;

				f64 factor = 2.0;
				if (k0 == k1)
					factor = 1.0;

				f64 tmn = pt2_cache[PT2_CACHE_INDEX(m-1,n-m)];
				const f64 delta_mn = (m == n) ? 1.0 : 0.0;

				f64 v_mn_pq = G0(&pt,component,component)*V_closed(&pt, hermite_cache, component,component, m,n,p,q);
				f64 tpq = pt2_cache[PT2_CACHE_INDEX(p-1,q-p)];

				const f64 delta_pq = (p == q) ? 1.0 : 0.0;
				const f64 coeff = 2.0 + c_root_2_minus_2*(delta_mn + delta_pq) + c_3_minus_2_root_2*(delta_mn*delta_pq);
				E_mn_pq += factor*coeff*tmn*tpq*v_mn_pq;
			}

		}

		sbmf_log_info("\t\tmn,pq: %.10e", E_mn_pq);

		E3 = E_00_00 + E_m0_n0 + E_mn_pq;
	}
	sbmf_log_info("\tE3: %e", E3);

	cudaFree(hermite_cache_device);

	return (struct pt_result) {
		.E0 = E0,
		.E1 = E1,
		.E2 = E2,
		.E3 = E3,
	};
}
