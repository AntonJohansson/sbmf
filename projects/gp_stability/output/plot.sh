#!/bin/sh

gnuplot -p -e "
	set term epslatex;
	set output 'imgs/gp_stability.eps';
	set size 0.75,0.75;
	set key at graph 0.95,0.2;
	set key box;
	set format y '%.1f';
	set grid;

	set xlabel 'Non-linear coupling constant \$\\lambda\$ (a.\,u.)';
	set ylabel 'Energy (a.\,u.)';
	set xrange [-10:5];
	plot
		'out_sgl' u 1:2 	w lines lw 2 title 'Harm. pot.',
		'out_dbl' u 1:2 	w lines lw 2 title 'Harm. pot. + Guassian'
"

gnuplot -p -e "
	set term epslatex;
	set output 'imgs/gp_stability_ex.eps';
	set size 0.75,0.75;
	set key at graph 0.975,0.975;
	set key box;
	set format y '%.1f';
	set grid;

	set xlabel '\$x\$ (a.\,u.)';
	set ylabel 'Density per particle (a.\,u.)';
	set xrange [-5:5];
	plot
		'ex_a' u 1:3 	w lines lw 2 title 'State 1',
		'ex_b' u 1:3 	w lines lw 2 title 'State 2'
"
