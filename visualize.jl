using DelimitedFiles, Plots
pyplot()

N = 128;

function load_wf()
	wf_data = readdlm("gss_item_wf");
	global wf = reshape(wf_data, N,N);

	#pt_data = readdlm("gss_item_pt");
	#global pt = reshape(pt_data, N,N);
end

function load_ev()
	global ev  = readdlm("gss_item_ev");
end

function load_vv()
	global vv = readdlm("gss_item_vv");
end

function load_grid()
	global grid = readdlm("debug_grid");
end

function load_all()
	load_wf();
	load_ev();
	load_vv();
end

function plot_wf()
	p1 = surface(wf);
end

function plot_vv(n)
	vmat = reshape(vv[n,:], N,N);
	surface(vmat);
end

function plot_grid()
	scatter(grid[:,1], grid[:,2]);
end

#err_data = readdlm("gss_item_error");
#plot(1:size(err_data, 1), err_data);
#
#dbg_data = readdlm("gss_item_debug");
#plot(1:size(err_data, 1), err_data);
