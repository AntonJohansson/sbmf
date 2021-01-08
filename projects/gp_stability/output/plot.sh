#!/bin/sh

	#set key box;
gnuplot -p -e "
	set size 0.75,0.75;
	set key at graph 0.95,0.5;
	set key spacing 1.5;
	set format y '%.1f';
	set term epslatex;
	set output 'imgs/gp_stability.eps';

	set xlabel 'Non-linear coupling constant \$\\lambda\$ (a.\,u.)';
	set ylabel 'Energy (a.\,u.)';
	set xrange [-10:10];
	plot
		'out_sgl_32' 		u 1:2 	w linespoints lw 2 title '\$V_\\mathrm{H}(x)\$',
		'out_sgl_32_hmix' 	u 1:2 	w linespoints lw 2 title '\$V_\\mathrm{H}(x)\$ w. \$\\op H\$ mix',
		'out_dbl_32' 		u 1:2 	w linespoints lw 2 title '\$V_\\mathrm{H,G}(x)\$',
		'out_dbl_32_hmix' 	u 1:2 	w linespoints lw 2 title '\$V_\\mathrm{H,G}(x)\$ w. \$\\op H\$ mix',
"

#gnuplot -p -e "
#	set term epslatex;
#	set output 'imgs/gp_stability_ex.eps';
#	set size 0.75,0.75;
#	set key at graph 0.975,0.975;
#	set key box;
#	set format y '%.1f';
#	set grid;
#
#	set xlabel '\$x\$ (a.\,u.)';
#	set ylabel 'Density per particle (a.\,u.)';
#	set xrange [-5:5];
#	plot
#		'ex_a' u 1:3 	w lines lw 2 title 'State 1',
#		'ex_b' u 1:3 	w lines lw 2 title 'State 2'
#"
