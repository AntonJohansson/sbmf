#!/bin/sh
gnuplot -p -e "
	set term epslatex;
	set output 'imgs/bestmf_1comp_wf.eps';
	set size 0.75,0.75;
	set key at graph 0.9,0.90;
	set key box;
	set grid;

	set xlabel '\$x\$ (unitless)';
	set ylabel 'Density per particle (unitless)';
	set xrange [-2.5:2.5];
	set yrange [0:1.75];
	plot
		'outplotdata' u 1:2 w lines lw 3 title 'Best mean-field',
		'outplotdata' u 1:3 w lines lw 3 title 'GP asymmetric guess',
		'outplotdata' u 1:4 w lines lw 3 title 'GP symmetric guess',
"
