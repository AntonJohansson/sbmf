__host__ __device__
static inline void map_to_triangular_index(u32 k, u32 N, u32* m, u32* n) {
	*m = k / N;
	*n = k % N;
	if (*m > *n) {
		*m = N - *m - 0;
		*n = N - *n - 1;
	}
}

__host__ __device__
static inline f64 V_closed(const f64* cache, const f64* phi_a, const f64* phi_b, const f64* phi_c, const f64* phi_d, const u32 size) {
	f64 sum = 0.0;
	for (u32 a = 0; a < size; ++a) {
		for (u32 b = 0; b < size; ++b) {
			for (u32 c = 0; c < size; ++c) {
				for (u32 d = 0; d < size; ++d) {
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

__global__
static void device_sum_reduction(f64* out, f64* arr, const u32 len) {
	f64 sum = 0.0;
	for (u32 i = 0; i < len; ++i) {
		sum += arr[i];
	}
	*out = sum;
}

enum pt_mode {
	MODE_RSPT = 0,
	MODE_ENPT = 1,
};

__global__
static void rspt_3_mnpq_1comp(enum pt_mode mode, f64 g, const u32 num_sb_states, const u32 num_mb_states, const u32 num_interactions, f64* pt2_cache, f64* hermite_cache, f64* coeffs, f64* output) {
	const f64 c_root_2_minus_2 = sqrt(2.0) - 2.0;
	const f64 c_3_minus_2_root_2 = 3.0 - 2.0*sqrt(2.0);

	const u32 k = blockIdx.x*blockDim.x + threadIdx.x;
	if (k >= num_interactions)
		return;
	u32 k0, k1;
	map_to_triangular_index(k, num_mb_states, &k0, &k1);

	f64 factor = 2.0;
	if (k0 == k1) {
		if (mode == MODE_RSPT)
			factor = 1.0;
		else if (mode == MODE_ENPT)
			return;
	}

	u32 m, n;
	map_to_triangular_index(k0, num_sb_states-1, &m, &n);
	m += 1; n += 1;

	u32 p, q;
	map_to_triangular_index(k1, num_sb_states-1, &p, &q);
	p += 1; q += 1;

	const f64 tmn = pt2_cache[index2_cuda(m-1,n-1)];
	const f64 tpq = pt2_cache[index2_cuda(p-1,q-1)];

	const f64 delta_mn = (m == n) ? 1.0 : 0.0;
	const f64 delta_pq = (p == q) ? 1.0 : 0.0;

	f64 v_mn_pq = g*V_closed(hermite_cache,
			&coeffs[m*num_sb_states],
			&coeffs[n*num_sb_states],
			&coeffs[p*num_sb_states],
			&coeffs[q*num_sb_states],
			num_sb_states);


	const f64 coeff = 2.0 + c_root_2_minus_2*(delta_mn + delta_pq) + c_3_minus_2_root_2*(delta_mn*delta_pq);
	output[k] = factor*coeff*tmn*tpq*v_mn_pq;
}

__global__
static void rspt_3_mnpq_2comp(enum pt_mode mode, f64 g, const u32 num_sb_states, const u32 num_mb_states, const u32 num_interactions, f64* pt2_cache, f64* hermite_cache, f64* coeffsA, f64* coeffsB, f64* output) {
	const u32 k = blockIdx.x*blockDim.x + threadIdx.x;
	if (k >= num_interactions)
		return;

	u32 k0, k1;
	map_to_triangular_index(k, num_mb_states, &k0, &k1);

	f64 factor = 2.0;
	if (k0 == k1) {
		if (mode == MODE_RSPT)
			factor = 1.0;
		else if (mode == MODE_ENPT) {
			output[k] = 0;
			return;
		}
	}

	const u32 m = k0 % (num_sb_states-1) + 1;
	const u32 n = k0 / (num_sb_states-1) + 1;
	const u32 p = k1 % (num_sb_states-1) + 1;
	const u32 q = k1 / (num_sb_states-1) + 1;

#define PT2_CACHE_INDEX(i, j) \
	(i)*(num_sb_states-1) + (j)

	const f64 tmn = pt2_cache[PT2_CACHE_INDEX(m-1, n-1)];
	const f64 tpq = pt2_cache[PT2_CACHE_INDEX(p-1, q-1)];
#undef PT2_CACHE_INDEX

	const f64 v_mn_pq = g * V_closed(hermite_cache,
			&coeffsA[m*num_sb_states],
			&coeffsB[n*num_sb_states],
			&coeffsA[p*num_sb_states],
			&coeffsB[q*num_sb_states],
			num_sb_states);
	output[k] = factor*v_mn_pq*tmn*tpq;
}

struct Vp_params {
	u32 coeff_count;
	f64* i;
	f64* j;
	nlse_operator_func* pert;
};

void Vp_integrand(f64* out, f64* in, u32 len, void* data);

static inline f64 en_nhn_new(f64* phi_m, f64* phi_n, const u32 num_sb_states, nlse_operator_func* pert) {
	f64 sum = 0;
	for (u32 k = 0; k < num_sb_states; ++k) {
		sum += phi_m[k]*phi_n[k]*ho_eigenval(k);

	}

	/*
	 * In the case that we're dealing with a perturbation to the
	 * basis potential, we need to compute <i|Vp|j> numerically
	 * with Vp being the pertubation
	 */

	if (pert) {
		struct Vp_params p = {
			.coeff_count = num_sb_states,
			.i = phi_m,
			.j = phi_n,
			.pert = pert,
		};

		struct quadgk_settings settings = {
			.gk = gk20,
			.abs_error_tol = 1e-15,
			.rel_error_tol = 1e-15,
			.max_iters = 500,
			.userdata = &p,
		};

		u8 quadgk_memory[quadgk_required_memory_size(&settings)];

		struct quadgk_result res;
		quadgk_infinite_interval(Vp_integrand, &settings, quadgk_memory, &res);
		assert(res.converged);

		sum += res.integral;
	}

	return sum;
}



/**************************************************************************************************************************************************/




static struct pt_result perturbation_theory_1comp(enum pt_mode mode, f64 g, i64 N, const f64* hermite_cache, const u32 hermite_cache_size, const struct eigen_result_real* states, const f64 groundstate_energy, const f64* double_subst_energy_diffs, const u32 num_sb_states) {
	f64* device_states;
	cudaMalloc(&device_states, num_sb_states*num_sb_states*sizeof(f64));
	cudaMemcpy(device_states, states->eigenvectors, num_sb_states*num_sb_states*sizeof(f64), cudaMemcpyHostToDevice);

	f64* hermite_cache_device;
	cudaMalloc(&hermite_cache_device, hermite_cache_size);
	cudaMemcpy(hermite_cache_device, hermite_cache, hermite_cache_size, cudaMemcpyHostToDevice);

	/* Zeroth order PT */
	sbmf_log_info("Starting zeroth order PT");
	f64 E0 = groundstate_energy;
	sbmf_log_info("\tE0: %e", E0);

	/* This particular integral shows up in zeroth and third order rspt */
	const f64 v_00_00 = V_closed(hermite_cache,
			&states->eigenvectors[0],
			&states->eigenvectors[0],
			&states->eigenvectors[0],
			&states->eigenvectors[0],
			num_sb_states);

	/* First order PT */
	f64 E1 = 0.0;
	if (mode == MODE_RSPT) {
		sbmf_log_info("Starting first order PT");
		E1 = -0.5*g*N*(N-1)*v_00_00;
		sbmf_log_info("\tE1: %e", E1);
	}

	const u32 pt2_cache_size = size2_cuda(num_sb_states-1);
	f64 pt2_cache[pt2_cache_size];

	/* Second order PT */
	sbmf_log_info("Starting second order PT");
	f64 E2 = 0.0;
	{
		/*
		 * Double substitutions (both excitations within same component),
		 * loop over unique pairs (j,k).
		 */
#pragma omp parallel for reduction(+: E2)
		for (u32 m = 1; m < num_sb_states; ++m) {
			for (u32 n = m; n < num_sb_states; ++n) {
				const f64 v_mn_00 = V_closed(hermite_cache,
						&states->eigenvectors[m*num_sb_states],
						&states->eigenvectors[n*num_sb_states],
						&states->eigenvectors[0*num_sb_states],
						&states->eigenvectors[0*num_sb_states],
						num_sb_states);
				const f64 factor = (m == n) ? 1.0/sqrt(2.0) : 1.0;
				const f64 me = factor*g*sqrt(N*(N-1))*v_mn_00;

				const f64 Ediff = double_subst_energy_diffs[index2_cuda(m-1,n-1)];

				pt2_cache[index2_cuda(m-1,n-1)] = me/Ediff;

				E2 += me*me/(Ediff);
			}
		}
	}
	sbmf_log_info("\tE2: %e", E2);

	f64* device_pt2_cache;
	cudaMalloc(&device_pt2_cache, pt2_cache_size*sizeof(f64));
	cudaMemcpy(device_pt2_cache, pt2_cache, pt2_cache_size*sizeof(f64), cudaMemcpyHostToDevice);

	/* Third order PT */
	sbmf_log_info("Starting third order PT");
	f64 E3 = 0.0;
	{
		f64 E_00_00 = 0;
		if (mode == MODE_RSPT) {
#pragma omp parallel for reduction(+: E_00_00)
			for (u32 m = 1; m < num_sb_states; ++m) {
				for (u32 n = m; n < num_sb_states; ++n) {
					const f64 tmn = pt2_cache[index2_cuda(m-1,n-1)];
					E_00_00 += tmn*tmn;
				}
			}

			E_00_00 *= g*v_00_00;
			sbmf_log_info("\t\t00,00: %.10e", E_00_00);
		}

		/* Number of many-body states, excludes 0,0 */
		const u32 num_mb_states = size2_cuda(num_sb_states-1);

		f64 E_m0_n0 = 0;
		{
			const f64 c_root_2_minus_1 = sqrt(2.0) - 1.0;
			const f64 c_3_minus_2_root_2 = 3.0 - 2.0*sqrt(2.0);

#pragma omp parallel for reduction(+: E_m0_n0)
			for (u32 k = 0; k < num_mb_states; ++k) {
				u32 m, n;
				map_to_triangular_index(k, num_sb_states-1, &m, &n);
				if (mode == MODE_ENPT && m == n)
					continue;

				m += 1;
				n += 1;

				const f64 v_m0_n0 = g*V_closed(hermite_cache,
						&states->eigenvectors[m*num_sb_states],
						&states->eigenvectors[0*num_sb_states],
						&states->eigenvectors[n*num_sb_states],
						&states->eigenvectors[0*num_sb_states],
						num_sb_states);

				f64 sum = 0.0;
				for (u32 p = 1; p < num_sb_states; ++p) {

					const f64 tmp = pt2_cache[index2_cuda(m-1,p-1)];
					const f64 tnp = pt2_cache[index2_cuda(n-1,p-1)];

					const f64 delta_mp = (m == p) ? 1.0 : 0.0;
					const f64 delta_np = (n == p) ? 1.0 : 0.0;

					const f64 coeff = 1 + c_root_2_minus_1*(delta_mp + delta_np) + c_3_minus_2_root_2*(delta_mp*delta_np);
					sum += coeff * tmp * tnp;
				}

				f64 factor = (m == n) ? 1.0 : 2.0;

				E_m0_n0 += factor * v_m0_n0 * sum;
			}

			E_m0_n0 *= (N - 3);
		}
		sbmf_log_info("\t\tm0,n0: %.10e", E_m0_n0);

		f64 E_mn_pq = 0;
		{
			const u32 num_interactions = (num_mb_states*(num_mb_states+1))/2;

			f64* device_output;
			cudaMalloc(&device_output, num_interactions*sizeof(f64));

			const u32 blocks = num_interactions/256 + 1;
			rspt_3_mnpq_1comp<<<blocks, 256>>>(
					mode,
					g,
					num_sb_states, num_mb_states, num_interactions,
					device_pt2_cache,
					hermite_cache_device,
					device_states,
					device_output
					);

			f64* res;
			cudaMalloc(&res, sizeof(f64));
			device_sum_reduction<<<1,1>>>(res, device_output, num_interactions);
			cudaMemcpy(&E_mn_pq, res, sizeof(f64), cudaMemcpyDeviceToHost);
			cudaFree(res);

			cudaFree(device_output);
		}

		sbmf_log_info("\t\tmn,pq: %.10e", E_mn_pq);

		E3 = E_00_00 + E_m0_n0 + E_mn_pq;
	}
	sbmf_log_info("\tE3: %e", E3);

	cudaFree(hermite_cache_device);
	cudaFree(device_states);
	cudaFree(device_pt2_cache);

	return (struct pt_result) {
		.E0 = E0,
		.E1 = E1,
		.E2 = E2,
		.E3 = E3,
	};
}

static struct pt_result perturbation_theory_2comp(enum pt_mode mode, f64 gAA, f64 gAB, i64 NA, i64 NB, const f64* hermite_cache, const u32 hermite_cache_size, const struct eigen_result_real* statesA, const struct eigen_result_real* statesB,
		const f64 groundstate_energy,
		const f64* double_subst_energy_diffs_AA,
		const f64* double_subst_energy_diffs_BB,
		const f64* double_subst_energy_diffs_AB,
		const u32 num_sb_states) {

	f64* device_states_A;
	f64* device_states_B;
	cudaMalloc(&device_states_A, num_sb_states*num_sb_states*sizeof(f64));
	cudaMemcpy(device_states_A, statesA->eigenvectors, num_sb_states*num_sb_states*sizeof(f64), cudaMemcpyHostToDevice);
	cudaMalloc(&device_states_B, num_sb_states*num_sb_states*sizeof(f64));
	cudaMemcpy(device_states_B, statesB->eigenvectors, num_sb_states*num_sb_states*sizeof(f64), cudaMemcpyHostToDevice);

	f64* hermite_cache_device;
	cudaMalloc(&hermite_cache_device, hermite_cache_size);
	cudaMemcpy(hermite_cache_device, hermite_cache, hermite_cache_size, cudaMemcpyHostToDevice);

	struct pt_result res_A = perturbation_theory_1comp(mode, gAA, NA, hermite_cache, hermite_cache_size, statesA, groundstate_energy, double_subst_energy_diffs_AA, num_sb_states);
	struct pt_result res_B = perturbation_theory_1comp(mode, gAA, NB, hermite_cache, hermite_cache_size, statesB, groundstate_energy, double_subst_energy_diffs_BB, num_sb_states);

	sbmf_log_info("Starting zeroth order PT");
	//f64 E0 = res_A.E0 + res_B.E0;
	f64 E0 = groundstate_energy;
	sbmf_log_info("\tE0: %e", E0);

	f64 E1 = 0.0;
	if (mode == MODE_RSPT) {
		sbmf_log_info("Starting first order PT");
		E1 = res_A.E1 + res_B.E1;
		E1 += -gAB*NA*NB*V_closed(hermite_cache,
				&statesA->eigenvectors[0],
				&statesB->eigenvectors[0],
				&statesA->eigenvectors[0],
				&statesB->eigenvectors[0],
				num_sb_states);
		sbmf_log_info("\tE1: %e", E1);
	}

	const u32 pt2_cache_size = (num_sb_states-1)*(num_sb_states-1);
	f64 pt2_cache[pt2_cache_size];

	/* Assumes i in [0,states_to_include), j in [0,states_to_include) */
#define PT2_CACHE_INDEX(i, j) \
	(i)*(num_sb_states-1) + (j)

	/* Second order PT */
	sbmf_log_info("Starting second order PT");
	f64 E2 = 0.0;
	{
		E2 = res_A.E2 + res_B.E2;
		/*
		 * Double substitutions (separate components).
		 * A,B refers to components.
		 */
#pragma omp parallel for reduction(+: E2)
		for (u32 m = 1; m < num_sb_states; ++m) {
			for (u32 n = 1; n < num_sb_states; ++n) {
				const f64 me = gAB*sqrt(NA*NB)
					* V_closed(hermite_cache,
							&statesA->eigenvectors[m*num_sb_states],
							&statesB->eigenvectors[n*num_sb_states],
							&statesA->eigenvectors[0*num_sb_states],
							&statesB->eigenvectors[0*num_sb_states],
							num_sb_states);

				const f64 Ediff = double_subst_energy_diffs_AB[(m-1)*(num_sb_states-1) + (n-1)];

				pt2_cache[PT2_CACHE_INDEX(m-1,n-1)] = me/Ediff;

				E2 += me*me/(Ediff);
			}
		}
	}
	sbmf_log_info("\tE2: %e", E2);

	f64* device_pt2_cache;
	cudaMalloc(&device_pt2_cache, pt2_cache_size*sizeof(f64));
	cudaMemcpy(device_pt2_cache, pt2_cache, pt2_cache_size*sizeof(f64), cudaMemcpyHostToDevice);

	sbmf_log_info("Starting third order PT");
	f64 E3 = 0.0;
	{
		E3 = res_A.E3 + res_B.E3;

		f64 E_00_00 = 0.0;
		if (mode == MODE_RSPT) {

			f64 sum = 0;
#pragma omp parallel for reduction(+: sum)
			for (u32 m = 1; m < num_sb_states; ++m) {
				for (u32 n = 1; n < num_sb_states; ++n) {
					const f64 tmn = pt2_cache[PT2_CACHE_INDEX(m-1, n-1)];
					sum += tmn*tmn;
				}
			}

			const f64 v_00_00 = V_closed(hermite_cache,
					&statesA->eigenvectors[0],
					&statesB->eigenvectors[0],
					&statesA->eigenvectors[0],
					&statesB->eigenvectors[0],
					num_sb_states);
			E_00_00 += gAB * v_00_00 * sum;
		}
		sbmf_log_info("\t\t00,00: %.10e", E_00_00);

		f64 E_m0_n0 = 0.0;
		{
			{
#pragma omp parallel for collapse(2) reduction(+: E_m0_n0)
				for (u32 m = 1; m < num_sb_states; ++m) {
					for (u32 n = 1; n < num_sb_states; ++n) {
						if (mode == MODE_ENPT && m == n)
							continue;

						f64 sumA = 0;
						f64 sumB = 0;
						for (u32 p = 1; p < num_sb_states; ++p) {
							const f64 tmp = pt2_cache[PT2_CACHE_INDEX(m-1, p-1)];
							const f64 tnp = pt2_cache[PT2_CACHE_INDEX(n-1, p-1)];
							sumA += tmp*tnp;

							const f64 tpm = pt2_cache[PT2_CACHE_INDEX(p-1, m-1)];
							const f64 tpn = pt2_cache[PT2_CACHE_INDEX(p-1, n-1)];
							sumB += tpm*tpn;
						}

						const f64 v_AA_m0_n0 = gAA*(NA-1)*V_closed(hermite_cache,
								&statesA->eigenvectors[m*num_sb_states],
								&statesA->eigenvectors[0*num_sb_states],
								&statesA->eigenvectors[n*num_sb_states],
								&statesA->eigenvectors[0*num_sb_states],
								num_sb_states);
						const f64 v_BB_m0_n0 = gAA*(NB-1)*V_closed(hermite_cache,
								&statesB->eigenvectors[m*num_sb_states],
								&statesB->eigenvectors[0*num_sb_states],
								&statesB->eigenvectors[n*num_sb_states],
								&statesB->eigenvectors[0*num_sb_states],
								num_sb_states);
						const f64 v_AB_m0_n0 = - gAB*V_closed(hermite_cache,
								&statesA->eigenvectors[m*num_sb_states],
								&statesB->eigenvectors[0*num_sb_states],
								&statesA->eigenvectors[n*num_sb_states],
								&statesB->eigenvectors[0*num_sb_states],
								num_sb_states);
						const f64 v_AB_0m_0n = - gAB*V_closed(hermite_cache,
								&statesA->eigenvectors[0*num_sb_states],
								&statesB->eigenvectors[m*num_sb_states],
								&statesA->eigenvectors[0*num_sb_states],
								&statesB->eigenvectors[n*num_sb_states],
								num_sb_states);
						E_m0_n0 += (v_AA_m0_n0 + v_AB_m0_n0) * sumA;
						E_m0_n0 += (v_BB_m0_n0 + v_AB_0m_0n) * sumB;
					}
				}
			}
		}
		sbmf_log_info("\t\tm0,n0: %.10e", E_m0_n0);

		f64 E_mn_pq = 0;
		{
			const u32 num_mb_states = (num_sb_states-1)*(num_sb_states-1);
			const u32 num_interactions = (num_mb_states*(num_mb_states+1))/2;

			f64* device_output;
			cudaMalloc(&device_output, num_interactions*sizeof(f64));

			const u32 blocks = num_interactions/256 + 1;
			rspt_3_mnpq_2comp<<<blocks, 256>>>(
					mode,
					gAB,
					num_sb_states, num_mb_states, num_interactions,
					device_pt2_cache,
					hermite_cache_device,
					device_states_A,
					device_states_B,
					device_output
					);

			f64* res;
			cudaMalloc(&res, sizeof(f64));
			device_sum_reduction<<<1,1>>>(res, device_output, num_interactions);
			cudaMemcpy(&E_mn_pq, res, sizeof(f64), cudaMemcpyDeviceToHost);
			cudaFree(res);

			cudaFree(device_output);
		}
		sbmf_log_info("\t\tmn,pq: %.10e", E_mn_pq);

		E3 += E_00_00 + E_m0_n0 + E_mn_pq;
	}
	sbmf_log_info("\tE3: %e", E3);

	cudaFree(device_pt2_cache);
	cudaFree(hermite_cache_device);
	cudaFree(device_states_A);
	cudaFree(device_states_B);

	return (struct pt_result) {
		.E0 = E0,
		.E1 = E1,
		.E2 = E2,
		.E3 = E3,
	};
}

/**************************************************************************************************************************************************/

struct pt_result rspt_1comp_cuda_new(struct nlse_settings* settings, struct nlse_result res, u32 component, f64 g, i64 N) {
	/*
	 * The number of single body (sb) states is equal to the number
	 * of coefficients which is equal to the basis size
	 */
	const u32 num_sb_states = res.coeff_count;

	/* Find all eigenstates and eigenenergies of the hamiltonian passed in */
	struct eigen_result_real states;
	states = find_eigenpairs_full_real(res.hamiltonian[component]);
	for (u32 j = 0; j < num_sb_states; ++j) {
		f64_normalize(&states.eigenvectors[j*num_sb_states], &states.eigenvectors[j*num_sb_states], num_sb_states);
	}

	const u64 hermite_integral_count = size4_cuda(num_sb_states);
	const u64 hermite_cache_size = sizeof(f64)*hermite_integral_count;
	u32 memory_marker = sbmf_stack_marker();
	f64* hermite_cache = (f64*)sbmf_stack_push(hermite_cache_size);
	{
		sbmf_log_info("Precomputing %ld hermite integrals", hermite_integral_count);
		for (u32 i = 0; i < num_sb_states; ++i) {
			for (u32 j = i; j < num_sb_states; ++j) {
				for (u32 k = j; k < num_sb_states; ++k) {
					for (u32 l = k; l < num_sb_states; ++l) {
						hermite_cache[index4_cuda(i,j,k,l)] = hermite_integral_4(i,j,k,l);
					}
				}
			}
		}
	}


	f64 groundstate_energy = N*states.eigenvalues[0];
	/* Energies of double substitution states including the zero states */
	f64 double_subst_energy_diffs[size2_cuda(num_sb_states-1)];
	for (u32 m = 1; m < num_sb_states; ++m) {
		for (u32 n = m; n < num_sb_states; ++n) {
			double_subst_energy_diffs[index2_cuda(m-1,n-1)] = 2*states.eigenvalues[0] - states.eigenvalues[m] - states.eigenvalues[n];
		}
	}

	struct pt_result ptres = perturbation_theory_1comp(MODE_RSPT, g, N, hermite_cache, hermite_cache_size, &states, groundstate_energy, double_subst_energy_diffs, num_sb_states);
	sbmf_stack_free_to_marker(memory_marker);

	return ptres;
}

struct pt_result enpt_1comp_cuda_new(struct nlse_settings* settings, struct nlse_result res, u32 component, f64 g, i64 N) {
	/*
	 * The number of single body (sb) states is equal to the number
	 * of coefficients which is equal to the basis size
	 */
	const u32 num_sb_states = res.coeff_count;

	/* Find all eigenstates and eigenenergies of the hamiltonian passed in */
	struct eigen_result_real states;
	states = find_eigenpairs_full_real(res.hamiltonian[component]);
	for (u32 j = 0; j < num_sb_states; ++j) {
		f64_normalize(&states.eigenvectors[j*num_sb_states], &states.eigenvectors[j*num_sb_states], num_sb_states);
	}

	const u64 hermite_integral_count = size4_cuda(num_sb_states);
	const u64 hermite_cache_size = sizeof(f64)*hermite_integral_count;
	u32 memory_marker = sbmf_stack_marker();
	f64* hermite_cache = (f64*)sbmf_stack_push(hermite_cache_size);
	{
		sbmf_log_info("Precomputing %ld hermite integrals", hermite_integral_count);
		for (u32 i = 0; i < num_sb_states; ++i) {
			for (u32 j = i; j < num_sb_states; ++j) {
				for (u32 k = j; k < num_sb_states; ++k) {
					for (u32 l = k; l < num_sb_states; ++l) {
						hermite_cache[index4_cuda(i,j,k,l)] = hermite_integral_4(i,j,k,l);
					}
				}
			}
		}
	}

	/* Energies of double substitution states including the zero states */
	f64 double_subst_energy_diffs[size2_cuda(num_sb_states-1)];

	const f64 v_00_00 = V_closed(hermite_cache,
			&states.eigenvectors[0],
			&states.eigenvectors[0],
			&states.eigenvectors[0],
			&states.eigenvectors[0],
			num_sb_states);
	f64 groundstate_energy = N*en_nhn_new(&states.eigenvectors[0*num_sb_states], &states.eigenvectors[0*num_sb_states], num_sb_states, settings->spatial_pot_perturbation) + 0.5*g*N*(N-1)*v_00_00;

	for (u32 m = 1; m < num_sb_states; ++m) {
		for (u32 n = m; n < num_sb_states; ++n) {
			const f64 v_mn_mn = V_closed(hermite_cache,
					&states.eigenvectors[m*num_sb_states],
					&states.eigenvectors[n*num_sb_states],
					&states.eigenvectors[m*num_sb_states],
					&states.eigenvectors[n*num_sb_states],
					num_sb_states);
			const f64 v_m0_m0 = V_closed(hermite_cache,
					&states.eigenvectors[m*num_sb_states],
					&states.eigenvectors[0*num_sb_states],
					&states.eigenvectors[m*num_sb_states],
					&states.eigenvectors[0*num_sb_states],
					num_sb_states);
			const f64 v_n0_n0 = (m == n) ? v_m0_m0 : V_closed(hermite_cache,
					&states.eigenvectors[n*num_sb_states],
					&states.eigenvectors[0*num_sb_states],
					&states.eigenvectors[n*num_sb_states],
					&states.eigenvectors[0*num_sb_states],
					num_sb_states);

			const f64 dmn = (m == n) ? 1.0 : 0.0;
			f64 energy =
				2*en_nhn_new(&states.eigenvectors[0*num_sb_states], &states.eigenvectors[0*num_sb_states], num_sb_states, settings->spatial_pot_perturbation)
				-en_nhn_new(&states.eigenvectors[m*num_sb_states], &states.eigenvectors[m*num_sb_states], num_sb_states, settings->spatial_pot_perturbation)
				-en_nhn_new(&states.eigenvectors[n*num_sb_states], &states.eigenvectors[n*num_sb_states], num_sb_states, settings->spatial_pot_perturbation)
				- g*((2.0-dmn)*v_mn_mn + 2.0*(N-2)*(v_m0_m0+v_n0_n0) - (2*N-3)*v_00_00);

			double_subst_energy_diffs[index2_cuda(m-1,n-1)] = energy;
		}
	}

	struct pt_result ptres = perturbation_theory_1comp(MODE_ENPT, g, N, hermite_cache, hermite_cache_size, &states, groundstate_energy, double_subst_energy_diffs, num_sb_states);
	sbmf_stack_free_to_marker(memory_marker);

	return ptres;
}

struct pt_result rspt_2comp_cuda_new(struct nlse_settings* settings, struct nlse_result res, u32 compA, u32 compB, f64 gAA, f64 gAB, i64 NA, i64 NB) {
	/*
	 * The number of single body (sb) states is equal to the number
	 * of coefficients which is equal to the basis size
	 */
	const u32 num_sb_states = res.coeff_count;

	/* Find all eigenstates and eigenenergies of the hamiltonian passed in */
	struct eigen_result_real statesA, statesB;
	statesA = find_eigenpairs_full_real(res.hamiltonian[compA]);
	statesB = find_eigenpairs_full_real(res.hamiltonian[compB]);
	for (u32 j = 0; j < num_sb_states; ++j) {
		f64_normalize(&statesA.eigenvectors[j*num_sb_states], &statesA.eigenvectors[j*num_sb_states], num_sb_states);
		f64_normalize(&statesB.eigenvectors[j*num_sb_states], &statesB.eigenvectors[j*num_sb_states], num_sb_states);
	}

	const u64 hermite_integral_count = size4_cuda(num_sb_states);
	const u64 hermite_cache_size = sizeof(f64)*hermite_integral_count;
	u32 memory_marker = sbmf_stack_marker();
	f64* hermite_cache = (f64*)sbmf_stack_push(hermite_cache_size);
	{
		sbmf_log_info("Precomputing %ld hermite integrals", hermite_integral_count);
		for (u32 i = 0; i < num_sb_states; ++i) {
			for (u32 j = i; j < num_sb_states; ++j) {
				for (u32 k = j; k < num_sb_states; ++k) {
					for (u32 l = k; l < num_sb_states; ++l) {
						hermite_cache[index4_cuda(i,j,k,l)] = hermite_integral_4(i,j,k,l);
					}
				}
			}
		}
	}

	const f64 groundstate_energy = NA*statesA.eigenvalues[0] + NB*statesB.eigenvalues[0];

	/* Energies of double substitution states including the zero states */
	f64 double_subst_energies_AA[size2_cuda(num_sb_states-1)];
	f64 double_subst_energies_BB[size2_cuda(num_sb_states-1)];
	for (u32 m = 1; m < num_sb_states; ++m) {
		for (u32 n = m; n < num_sb_states; ++n) {
			double_subst_energies_AA[index2_cuda(m-1,n-1)] = 2*statesA.eigenvalues[0] - statesA.eigenvalues[m] - statesA.eigenvalues[n];
			double_subst_energies_BB[index2_cuda(m-1,n-1)] = 2*statesB.eigenvalues[0] - statesB.eigenvalues[m] - statesB.eigenvalues[n];
		}
	}

	f64 double_subst_energies_AB[(num_sb_states-1)*(num_sb_states-1)];
	for (u32 m = 1; m < num_sb_states; ++m) {
		for (u32 n = 1; n < num_sb_states; ++n) {
			double_subst_energies_AB[(m-1)*(num_sb_states-1) + (n-1)] = statesA.eigenvalues[0] + statesB.eigenvalues[0] - statesA.eigenvalues[m] - statesB.eigenvalues[n];
		}
	}

	struct pt_result ptres = perturbation_theory_2comp(MODE_RSPT, gAA, gAB, NA, NB, hermite_cache, hermite_cache_size, &statesA, &statesB, groundstate_energy, double_subst_energies_AA, double_subst_energies_BB, double_subst_energies_AB, num_sb_states);
	sbmf_stack_free_to_marker(memory_marker);

	return ptres;
}

struct pt_result enpt_2comp_cuda_new(struct nlse_settings* settings, struct nlse_result res, u32 compA, u32 compB, f64 gAA, f64 gAB, i64 NA, i64 NB) {
	/*
	 * The number of single body (sb) states is equal to the number
	 * of coefficients which is equal to the basis size
	 */
	const u32 num_sb_states = res.coeff_count;

	/* Find all eigenstates and eigenenergies of the hamiltonian passed in */
	struct eigen_result_real statesA, statesB;
	statesA = find_eigenpairs_full_real(res.hamiltonian[compA]);
	statesB = find_eigenpairs_full_real(res.hamiltonian[compB]);
	for (u32 j = 0; j < num_sb_states; ++j) {
		f64_normalize(&statesA.eigenvectors[j*num_sb_states], &statesA.eigenvectors[j*num_sb_states], num_sb_states);
		f64_normalize(&statesB.eigenvectors[j*num_sb_states], &statesB.eigenvectors[j*num_sb_states], num_sb_states);
	}

	const u64 hermite_integral_count = size4_cuda(num_sb_states);
	const u64 hermite_cache_size = sizeof(f64)*hermite_integral_count;
	u32 memory_marker = sbmf_stack_marker();
	f64* hermite_cache = (f64*)sbmf_stack_push(hermite_cache_size);
	{
		sbmf_log_info("Precomputing %ld hermite integrals", hermite_integral_count);
		for (u32 i = 0; i < num_sb_states; ++i) {
			for (u32 j = i; j < num_sb_states; ++j) {
				for (u32 k = j; k < num_sb_states; ++k) {
					for (u32 l = k; l < num_sb_states; ++l) {
						hermite_cache[index4_cuda(i,j,k,l)] = hermite_integral_4(i,j,k,l);
					}
				}
			}
		}
	}

	const f64 v_AA_00_00 = V_closed(hermite_cache,
			&statesA.eigenvectors[0],
			&statesA.eigenvectors[0],
			&statesA.eigenvectors[0],
			&statesA.eigenvectors[0],
			num_sb_states);
	const f64 v_BB_00_00 = V_closed(hermite_cache,
			&statesB.eigenvectors[0],
			&statesB.eigenvectors[0],
			&statesB.eigenvectors[0],
			&statesB.eigenvectors[0],
			num_sb_states);
	const f64 v_AB_00_00 = V_closed(hermite_cache,
			&statesA.eigenvectors[0],
			&statesB.eigenvectors[0],
			&statesA.eigenvectors[0],
			&statesB.eigenvectors[0],
			num_sb_states);

	const f64 groundstate_energy =
		  NA*en_nhn_new(&statesA.eigenvectors[0], &statesA.eigenvectors[0], num_sb_states, settings->spatial_pot_perturbation) + 0.5*gAA*NA*(NA-1)*v_AA_00_00
		+ NB*en_nhn_new(&statesB.eigenvectors[0], &statesB.eigenvectors[0], num_sb_states, settings->spatial_pot_perturbation) + 0.5*gAA*NB*(NB-1)*v_BB_00_00
		+ gAB*NA*NB*v_AB_00_00;

	/* Energies of double substitution states including the zero states */
	f64 double_subst_energy_diffs_AA[size2_cuda(num_sb_states-1)];
	f64 double_subst_energy_diffs_BB[size2_cuda(num_sb_states-1)];
	for (u32 m = 1; m < num_sb_states; ++m) {
		for (u32 n = m; n < num_sb_states; ++n) {
			f64 delta_mn = (m == n) ? 1.0 : 0.0;
			const f64 v_AA_mn_mn = V_closed(hermite_cache,
					&statesA.eigenvectors[m*num_sb_states],
					&statesA.eigenvectors[n*num_sb_states],
					&statesA.eigenvectors[m*num_sb_states],
					&statesA.eigenvectors[n*num_sb_states],
					num_sb_states);
			const f64 v_AA_m0_m0 = V_closed(hermite_cache,
					&statesA.eigenvectors[m*num_sb_states],
					&statesA.eigenvectors[0*num_sb_states],
					&statesA.eigenvectors[m*num_sb_states],
					&statesA.eigenvectors[0*num_sb_states],
					num_sb_states);
			const f64 v_AA_n0_n0 = (m == n) ? v_AA_m0_m0 : V_closed(hermite_cache,
					&statesA.eigenvectors[n*num_sb_states],
					&statesA.eigenvectors[0*num_sb_states],
					&statesA.eigenvectors[n*num_sb_states],
					&statesA.eigenvectors[0*num_sb_states],
					num_sb_states);
			const f64 v_AB_m0_m0 = V_closed(hermite_cache,
					&statesA.eigenvectors[m*num_sb_states],
					&statesB.eigenvectors[0*num_sb_states],
					&statesA.eigenvectors[m*num_sb_states],
					&statesB.eigenvectors[0*num_sb_states],
					num_sb_states);
			const f64 v_AB_n0_n0 = (m == n) ? v_AB_m0_m0 : V_closed(hermite_cache,
					&statesA.eigenvectors[n*num_sb_states],
					&statesB.eigenvectors[0*num_sb_states],
					&statesA.eigenvectors[n*num_sb_states],
					&statesB.eigenvectors[0*num_sb_states],
					num_sb_states);
			double_subst_energy_diffs_AA[index2_cuda(m-1,n-1)] =
				  2*en_nhn_new(&statesA.eigenvectors[0], &statesA.eigenvectors[0], num_sb_states, settings->spatial_pot_perturbation)
				  - en_nhn_new(&statesA.eigenvectors[m*num_sb_states], &statesA.eigenvectors[m*num_sb_states], num_sb_states, settings->spatial_pot_perturbation)
				  - en_nhn_new(&statesA.eigenvectors[n*num_sb_states], &statesA.eigenvectors[n*num_sb_states], num_sb_states, settings->spatial_pot_perturbation)
				  - gAA*((2.0-delta_mn)*v_AA_mn_mn + 2.0*(NA-2)*(v_AA_m0_m0 + v_AA_n0_n0) - (2.0*NA-3.0)*v_AA_00_00)
				  - gAB*NB*(v_AB_m0_m0 + v_AB_n0_n0 - 2.0*v_AB_00_00);


			const f64 v_BB_mn_mn = V_closed(hermite_cache,
					&statesB.eigenvectors[m*num_sb_states],
					&statesB.eigenvectors[n*num_sb_states],
					&statesB.eigenvectors[m*num_sb_states],
					&statesB.eigenvectors[n*num_sb_states],
					num_sb_states);
			const f64 v_BB_m0_m0 = V_closed(hermite_cache,
					&statesB.eigenvectors[m*num_sb_states],
					&statesB.eigenvectors[0*num_sb_states],
					&statesB.eigenvectors[m*num_sb_states],
					&statesB.eigenvectors[0*num_sb_states],
					num_sb_states);
			const f64 v_BB_n0_n0 = (m == n) ? v_AA_m0_m0 : V_closed(hermite_cache,
					&statesB.eigenvectors[n*num_sb_states],
					&statesB.eigenvectors[0*num_sb_states],
					&statesB.eigenvectors[n*num_sb_states],
					&statesB.eigenvectors[0*num_sb_states],
					num_sb_states);
			const f64 v_AB_0m_0m = V_closed(hermite_cache,
					&statesA.eigenvectors[0*num_sb_states],
					&statesB.eigenvectors[m*num_sb_states],
					&statesA.eigenvectors[0*num_sb_states],
					&statesB.eigenvectors[m*num_sb_states],
					num_sb_states);
			const f64 v_AB_0n_0n = (m == n) ? v_AB_m0_m0 : V_closed(hermite_cache,
					&statesA.eigenvectors[0*num_sb_states],
					&statesB.eigenvectors[n*num_sb_states],
					&statesA.eigenvectors[0*num_sb_states],
					&statesB.eigenvectors[n*num_sb_states],
					num_sb_states);
			double_subst_energy_diffs_BB[index2_cuda(m-1,n-1)] =
				  2*en_nhn_new(&statesB.eigenvectors[0], &statesB.eigenvectors[0], num_sb_states, settings->spatial_pot_perturbation)
				  - en_nhn_new(&statesB.eigenvectors[m*num_sb_states], &statesB.eigenvectors[m*num_sb_states], num_sb_states, settings->spatial_pot_perturbation)
				  - en_nhn_new(&statesB.eigenvectors[n*num_sb_states], &statesB.eigenvectors[n*num_sb_states], num_sb_states, settings->spatial_pot_perturbation)
				  - gAA*((2.0-delta_mn)*v_BB_mn_mn + 2.0*(NB-2)*(v_BB_m0_m0 + v_BB_n0_n0) - (2.0*NB-3.0)*v_BB_00_00)
				  - gAB*NB*(v_AB_0m_0m + v_AB_0n_0n - 2.0*v_AB_00_00);
		}
	}

	f64 double_subst_energy_diffs_AB[(num_sb_states-1)*(num_sb_states-1)];
	for (u32 m = 1; m < num_sb_states; ++m) {
		for (u32 n = 1; n < num_sb_states; ++n) {
			const f64 v_AA_m0_m0 = V_closed(hermite_cache,
					&statesA.eigenvectors[m*num_sb_states],
					&statesA.eigenvectors[0*num_sb_states],
					&statesA.eigenvectors[m*num_sb_states],
					&statesA.eigenvectors[0*num_sb_states],
					num_sb_states);
			const f64 v_BB_n0_n0 = V_closed(hermite_cache,
					&statesB.eigenvectors[n*num_sb_states],
					&statesB.eigenvectors[0*num_sb_states],
					&statesB.eigenvectors[n*num_sb_states],
					&statesB.eigenvectors[0*num_sb_states],
					num_sb_states);
			const f64 v_AB_mn_mn = V_closed(hermite_cache,
					&statesA.eigenvectors[m*num_sb_states],
					&statesB.eigenvectors[n*num_sb_states],
					&statesA.eigenvectors[m*num_sb_states],
					&statesB.eigenvectors[n*num_sb_states],
					num_sb_states);
			const f64 v_AB_m0_m0 = V_closed(hermite_cache,
					&statesA.eigenvectors[m*num_sb_states],
					&statesB.eigenvectors[0*num_sb_states],
					&statesA.eigenvectors[m*num_sb_states],
					&statesB.eigenvectors[0*num_sb_states],
					num_sb_states);
			const f64 v_AB_0n_0n = V_closed(hermite_cache,
					&statesA.eigenvectors[0*num_sb_states],
					&statesB.eigenvectors[n*num_sb_states],
					&statesA.eigenvectors[0*num_sb_states],
					&statesB.eigenvectors[n*num_sb_states],
					num_sb_states);
			double_subst_energy_diffs_AB[(m-1)*(num_sb_states-1) + (n-1)] =
				    en_nhn_new(&statesA.eigenvectors[0], &statesA.eigenvectors[0], num_sb_states, settings->spatial_pot_perturbation)
				  - en_nhn_new(&statesA.eigenvectors[m*num_sb_states], &statesA.eigenvectors[m*num_sb_states], num_sb_states, settings->spatial_pot_perturbation)
				  - gAA*(2.0*(NA-1)*v_AA_m0_m0 - (NA-1)*v_AA_00_00)
				  + en_nhn_new(&statesB.eigenvectors[0], &statesB.eigenvectors[0], num_sb_states, settings->spatial_pot_perturbation)
				  - en_nhn_new(&statesB.eigenvectors[n*num_sb_states], &statesB.eigenvectors[n*num_sb_states], num_sb_states, settings->spatial_pot_perturbation)
				  - gAA*(2.0*(NB-1)*v_BB_n0_n0 - (NB-1)*v_BB_00_00)
				  - gAB*(v_AB_mn_mn + (NB-1)*v_AB_m0_m0 + (NA-1)*v_AB_0n_0n - (NA + NB - 1)*v_AB_00_00);
		}
	}

	struct pt_result ptres = perturbation_theory_2comp(MODE_ENPT, gAA, gAB, NA, NB, hermite_cache, hermite_cache_size, &statesA, &statesB, groundstate_energy, double_subst_energy_diffs_AA, double_subst_energy_diffs_BB, double_subst_energy_diffs_AB, num_sb_states);
	sbmf_stack_free_to_marker(memory_marker);

	return ptres;
}
