typedef f64 coord_transform_input_func(f64 x, f64 a, f64 b);
typedef f64 coord_transform_output_func(f64 x);

struct gk_data gk7 = {
	.kronod_nodes = {
		9.9145537112081263920685469752598e-01,
		9.4910791234275852452618968404809e-01,
		8.6486442335976907278971278864098e-01,
		7.415311855993944398638647732811e-01,
		5.8608723546769113029414483825842e-01,
		4.0584515137739716690660641207707e-01,
		2.0778495500789846760068940377309e-01,
		0.0
	},
	.kronod_weights = {
		2.2935322010529224963732008059913e-02,
		6.3092092629978553290700663189093e-02,
		1.0479001032225018383987632254189e-01,
		1.4065325971552591874518959051021e-01,
		1.6900472663926790282658342659795e-01,
		1.9035057806478540991325640242055e-01,
		2.0443294007529889241416199923466e-01,
		2.0948214108472782801299917489173e-01
	},
	.gauss_weights = {
		1.2948496616886969327061143267787e-01,
		2.797053914892766679014677714229e-01,
		3.8183005050511894495036977548818e-01,
		4.1795918367346938775510204081658e-01
	},
	.kronod_size = 8,
	.gauss_size = 4,
};

struct gk_data gk10 = {
	.kronod_nodes = {
		9.956571630258080807355272806890028e-01,
		9.739065285171717200779640120844521e-01,
		9.301574913557082260012071800595083e-01,
		8.650633666889845107320966884234930e-01,
		7.808177265864168970637175783450424e-01,
		6.794095682990244062343273651148736e-01,
		5.627571346686046833390000992726941e-01,
		4.333953941292471907992659431657842e-01,
		2.943928627014601981311266031038656e-01,
		1.488743389816312108848260011297200e-01,
		0.0,
	},
	.kronod_weights = {
		1.169463886737187427806439606219205e-02,
		3.255816230796472747881897245938976e-02,
		5.475589657435199603138130024458018e-02,
		7.503967481091995276704314091619001e-02,
		9.312545458369760553506546508336634e-02,
		1.093871588022976418992105903258050e-01,
		1.234919762620658510779581098310742e-01,
		1.347092173114733259280540017717068e-01,
		1.427759385770600807970942731387171e-01,
		1.477391049013384913748415159720680e-01,
		1.494455540029169056649364683898212e-01,
	},
	.gauss_weights = {
		6.667134430868813759356880989333179e-02,
		1.494513491505805931457763396576973e-01,
		2.190863625159820439955349342281632e-01,
		2.692667193099963550912269215694694e-01,
		2.955242247147528701738929946513383e-01,
	},
	.kronod_size = 11,
	.gauss_size = 5,
};

struct gk_data gk15 = {
	.kronod_nodes = {
		9.980022986933970602851728401522712e-01,
		9.879925180204854284895657185866126e-01,
		9.677390756791391342573479787843372e-01,
		9.372733924007059043077589477102095e-01,
		8.972645323440819008825096564544959e-01,
		8.482065834104272162006483207742169e-01,
		7.904185014424659329676492948179473e-01,
		7.244177313601700474161860546139380e-01,
		6.509967412974169705337358953132747e-01,
		5.709721726085388475372267372539106e-01,
		4.850818636402396806936557402323506e-01,
		3.941513470775633698972073709810455e-01,
		2.991800071531688121667800242663890e-01,
		2.011940939974345223006283033945962e-01,
		1.011420669187174990270742314473923e-01,
		0.000000000000000000000000000000000e+00,
	},
	.kronod_weights = {
		5.377479872923348987792051430127650e-03,
		1.500794732931612253837476307580727e-02,
		2.546084732671532018687400101965336e-02,
		3.534636079137584622203794847836005e-02,
		4.458975132476487660822729937327969e-02,
		5.348152469092808726534314723943030e-02,
		6.200956780067064028513923096080293e-02,
		6.985412131872825870952007709914748e-02,
		7.684968075772037889443277748265901e-02,
		8.308050282313302103828924728610379e-02,
		8.856444305621177064727544369377430e-02,
		9.312659817082532122548687274734572e-02,
		9.664272698362367850517990762758934e-02,
		9.917359872179195933239317348460313e-02,
		1.007698455238755950449466626175697e-01,
		1.013300070147915490173747927674925e-01,
	},
	.gauss_weights = {
    	3.075324199611726835462839357720442e-02,
        7.036604748810812470926741645066734e-02,
        1.071592204671719350118695466858693e-01,
        1.395706779261543144478047945110283e-01,
        1.662692058169939335532008604812088e-01,
        1.861610000155622110268005618664228e-01,
        1.984314853271115764561183264438393e-01,
 		2.025782419255612728806201999675193e-01,
	},
	.kronod_size = 16,
	.gauss_size = 8,
};

struct gk_data gk20 = {
	.kronod_nodes = {
		9.988590315882776638383155765458630e-01,
		9.931285991850949247861223884713203e-01,
		9.815078774502502591933429947202169e-01,
		9.639719272779137912676661311972772e-01,
		9.408226338317547535199827222124434e-01,
		9.122344282513259058677524412032981e-01,
		8.782768112522819760774429951130785e-01,
		8.391169718222188233945290617015207e-01,
		7.950414288375511983506388332727879e-01,
		7.463319064601507926143050703556416e-01,
		6.932376563347513848054907118459315e-01,
		6.360536807265150254528366962262859e-01,
		5.751404468197103153429460365864251e-01,
		5.108670019508270980043640509552510e-01,
		4.435931752387251031999922134926401e-01,
		3.737060887154195606725481770249272e-01,
		3.016278681149130043205553568585923e-01,
		2.277858511416450780804961953685746e-01,
		1.526054652409226755052202410226775e-01,
		7.652652113349733375464040939883821e-02,
		0.000000000000000000000000000000000e+00
	},
	.kronod_weights = {
		3.073583718520531501218293246030987e-03,
		8.600269855642942198661787950102347e-03,
		1.462616925697125298378796030886836e-02,
		2.038837346126652359801023143275471e-02,
		2.588213360495115883450506709615314e-02,
		3.128730677703279895854311932380074e-02,
		3.660016975820079803055724070721101e-02,
		4.166887332797368626378830593689474e-02,
		4.643482186749767472023188092610752e-02,
		5.094457392372869193270767005034495e-02,
		5.519510534828599474483237241977733e-02,
		5.911140088063957237496722064859422e-02,
		6.265323755478116802587012217425498e-02,
		6.583459713361842211156355696939794e-02,
		6.864867292852161934562341188536780e-02,
		7.105442355344406830579036172321017e-02,
		7.303069033278666749518941765891311e-02,
		7.458287540049918898658141836248753e-02,
		7.570449768455667465954277537661656e-02,
		7.637786767208073670550283503806100e-02,
		7.660071191799965644504990153010174e-02
	},
	.gauss_weights = {
		1.761400713915211831186196235185282e-02,
		4.060142980038694133103995227493211e-02,
		6.267204833410906356950653518704161e-02,
		8.327674157670474872475814322204621e-02,
		1.019301198172404350367501354803499e-01,
		1.181945319615184173123773777113823e-01,
		1.316886384491766268984944997481631e-01,
		1.420961093183820513292983250671649e-01,
		1.491729864726037467878287370019694e-01,
		1.527533871307258506980843319550976e-01
	},
	.kronod_size = 21,
	.gauss_size = 10,
};






























struct segment {
	f64 integral;
	f64 error;

	f64 start;
	f64 end;
};

static inline f64 norm(f64 r) {
	return r*r; //fabs(r);
}

struct eval_result {
	struct segment seg;
	u32 func_evals;
	bool valid;
};

static void evaluate_rule(integrand_func* f, f64 start, f64 end, const struct quadgk_settings* settings, struct eval_result* res) {
	const struct gk_data* gk = &settings->gk;
	f64 mid = 0.5 * (end - start);

	/*
	 *  	0 -- if even size
	 * 	1 -- if odd size
	 */
	u32 is_order_odd = 1 - (gk->kronod_size & 1);

	/*
	 * fg/k - gauss/kronod function eval
	 * Ig/k - gauss/kronod integral estimate
	 * xg/k - gauss/kronod sample nodes
	 */

	/* Size of the gauss weights array (other than the last element if odd order) */
	u32 size = gk->gauss_size - is_order_odd;

	/* Compute all sample points */
	u32 sample_size = 4*size + 3;
	f64 sample_points[sample_size];
	f64 transformed_sample_points[sample_size];
	for (u32 i = 0; i < size; ++i) {
		f64 xg = gk->kronod_nodes[2*i+1];
		f64 xk = gk->kronod_nodes[2*i];

		sample_points[4*i + 0] = start + (1-xg)*mid;
		sample_points[4*i + 1] = start + (1+xg)*mid;
		sample_points[4*i + 2] = start + (1-xk)*mid;
		sample_points[4*i + 3] = start + (1+xk)*mid;
	}
	sample_points[4*size + 0] = start + mid;
	sample_points[4*size + 1] = start + (1-gk->kronod_nodes[gk->kronod_size-2])*mid;
	sample_points[4*size + 2] = start + (1+gk->kronod_nodes[gk->kronod_size-2])*mid;

	/* Apply transform scaling to input */
	for (u32 i = 0; i < sample_size; ++i) {
		transformed_sample_points[i] = sample_points[i]/(1.0 - sample_points[i]*sample_points[i]);
	}

	/* Sample the function */
	f64 output[sample_size];
	f(output, transformed_sample_points, sample_size, settings->userdata);

	/* Apply transform output scaling */
	f64 one_minus_t2, scale;
	for (u32 i = 0; i < sample_size; ++i) {
		one_minus_t2 = (1.0 - sample_points[i]*sample_points[i]);
		scale = (1.0 + sample_points[i]*sample_points[i]) / (one_minus_t2*one_minus_t2);
		output[i] *= scale;
	}

	/* Compute the integral estimate */
	f64 Ig = 0, Ik = 0;
	for (u32 i = 0; i < size; ++i) {
		f64 fg = output[4*i + 0] + output[4*i + 1];
		f64 fk = output[4*i + 2] + output[4*i + 3];

		Ig += gk->gauss_weights[i] * fg;
		Ik += gk->kronod_weights[2*i+1] * fg + gk->kronod_weights[2*i] * fk;
	}

	// In the odd-order case, the last point has to be handled
	// differently
	if (is_order_odd  == 0) {
		Ik += gk->kronod_weights[gk->kronod_size-1] * output[4*size + 0];
	} else {
		Ig += gk->gauss_weights[gk->gauss_size-1] * output[4*size + 0];
		Ik += gk->kronod_weights[gk->kronod_size-1] * output[4*size + 0]
			+ gk->kronod_weights[gk->kronod_size-2] * (output[4*size + 1] + output[4*size + 2]);
	}

	f64 Ik_s = Ik*mid;
	f64 Ig_s = Ig*mid;
	f64 error = norm(Ik_s - Ig_s);

	bool valid = (!isinf(error) && !isnan(error));
	if (!valid) {
		sbmf_log_error("evaluation of gk rule resulted in [%s]", (isinf(error) == 1) ? "inf" : "nan");
	}

	/* Setup return values */
	res->seg = (struct segment){Ik_s, error, start, end};
	res->func_evals = 4*gk->gauss_size + 1;
	res->valid = valid;
}

/* Compare function for priority queue, sort segments by error in descending order. */
static bool compare_segments(void* a, void* b) {
	return (((struct segment*)a)->error > ((struct segment*)b)->error);
}

static inline bool should_exit(const struct quadgk_settings* settings, struct quadgk_result* res) {
	return (res->error <= settings->abs_error_tol ||
			res->error <= settings->rel_error_tol*norm(res->integral) ||
			res->performed_evals >= settings->max_evals);
}

void quadgk_infinite_interval(integrand_func* f, const struct quadgk_settings* settings, struct quadgk_result* res) {
	if (!settings) {
		sbmf_log_error("You need to provide a quadgk_settings struct!");
		return;
	}

	if (!res) {
		sbmf_log_error("You need to provide a quadgk_result struct!");
		return;
	}

	if (settings->gk.gauss_size == 0) {
		sbmf_log_error("You need to specify a gk rule!");
		return;
	}

	/* Interval for for transformed integrand */
	f64 start = -1.0, end = 1.0;

	res->integral = 0;
	res->error = 0;
	res->performed_evals = 0;
	res->converged = false;

	/*
	 * Perform one evaluation of the rule and see if we can
	 * exit early without having to initalize extra memory and
	 * so on.
	 */
	struct eval_result eval_res = {0};
	evaluate_rule(f, start, end, settings, &eval_res);
	res->performed_evals += eval_res.func_evals;

	struct segment s = eval_res.seg;
	res->integral = s.integral;
	res->error = s.error;

	if (!eval_res.valid || should_exit(settings, res)) {
		if (eval_res.valid && res->performed_evals <= settings->max_evals)
			res->converged = true;
		return;
	}

	/*
	 * If we get to this point we were not able to exit early and have do to more subdivisions
	 * to get desired error tolerance.
	 */

	u32 memory_marker = sbmf_stack_marker();

	struct prioqueue* pq = prioqueue_new(64, sizeof(struct segment), compare_segments);
	prioqueue_push(pq, &s);

	struct segment largest_error_seg = {0};
	struct eval_result left_eval_res = {0}, right_eval_res = {0};

	while (!should_exit(settings, res)) {
		prioqueue_pop(pq, &largest_error_seg);

		f64 midpoint = 0.5 * (largest_error_seg.start + largest_error_seg.end);
		evaluate_rule(f, largest_error_seg.start, midpoint, settings, &left_eval_res);
		evaluate_rule(f, midpoint, largest_error_seg.end, settings, &right_eval_res);
		if (!left_eval_res.valid || !right_eval_res.valid)
			return;

		res->performed_evals += left_eval_res.func_evals + right_eval_res.func_evals;
		res->integral = (res->integral - largest_error_seg.integral) + left_eval_res.seg.integral + right_eval_res.seg.integral;
		res->error = (res->error - largest_error_seg.error) + left_eval_res.seg.error + right_eval_res.seg.error;

		prioqueue_push(pq, &left_eval_res.seg);
		prioqueue_push(pq, &right_eval_res.seg);
	}

	sbmf_stack_free_to_marker(memory_marker);

	if (res->performed_evals <= settings->max_evals)
		res->converged = true;
}
