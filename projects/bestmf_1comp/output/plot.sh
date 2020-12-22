#!/bin/sh

gnuplot -p -e "
	set term epslatex;
	set output 'imgs/bestmf_1comp_wf.eps';
	set size 0.75,0.75;
	set key at graph 0.9,0.90;
	set key box;
	set grid;
	set format y '%.1f';

	set xlabel '\$x\$ (a.u.)';
	set ylabel 'Density per particle (unitless)';
	set xrange [-2.5:2.5];
	set yrange [0:1.75];
	plot
		'outplotdata' u 1:2 w lines lw 4 title 'Best mean-field',
		'outplotdata' u 1:3 w lines lw 4 title 'GP symmetric guess',
		'outplotdata' u 1:4 w lines lw 4 lt 4 title 'GP asymmetric guess',
"
