#!/bin/sh

gnuplot -p -e "
	set term epslatex;
	set output 'imgs/gp_pt_cisd.eps';
	set size 0.75,0.75;
	set key at graph 0.9,0.90;
	set key spacing 1.5;
	set format y '%.3f';

	set xlabel 'Basis size (a.\,u.)';
	set ylabel 'Energy (a.\,u.)';
	plot
		'out' u log(1):2:xtic(1) 	w linespoints lw 2 title '\$E_\\mathrm{GP}\$',
		'out' u log(1):7 			w linespoints lw 2 title '\$E_\\mathrm{RSPT}\$',
		'out_en' u log(1):7 			w linespoints lw 2 title '\$E_\\mathrm{ENPT}\$',
		2.713319038091934 lw 3 title '\$E_\\mathrm{CI}\$'
"
