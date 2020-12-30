#!/bin/sh

gnuplot -p -e "
	set term epslatex;
	set output 'imgs/gp_pt_cisd.eps';
	set size 0.75,0.75;
	set key at graph 0.9,0.90;
	set key box;
	set format y '%.3f';
	set grid;

	set xlabel 'Basis size (unitless)';
	set ylabel 'Energy (a.u.)';
	plot
		'out' u log(1):2:xtic(1) 	w linespoints lw 2 title '\$E_\\mathrm{GP}\$',
		'out' u log(1):7 			w linespoints lw 2 title '\$E_\\mathrm{PT}\$',
		2.71298958 lt 4 lw 3 title '\$E_\\mathrm{CI}\$'
"
