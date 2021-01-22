#!/bin/sh

gnuplot -p -e "
	set term epslatex;
	set output 'imgs/gp_pt_cisd.eps';
	set size 0.75,0.75;
	set key at graph 1.0,0.975;
	set key spacing 1.5;
	set format y '%.3f';

	set xlabel 'Basis size (a.\,u.)';
	set ylabel 'Energy (a.\,u.)';
	plot
		'out_rs' u log(1):2:xtic(1) 	w linespoints lw 2 title '\$E_\\mathrm{GP}\$',
		'out_rs' u log(1):3 			w linespoints dt 2 lc 4 lw 2 title '\$E_\\mathrm{RSPT2}\$',
		'out_rs' u log(1):4 			w linespoints lc 4 lw 2 title '\$E_\\mathrm{RSPT3}\$',
		'out_en' u log(1):3 			w linespoints dt 2 lc 6 lw 2 title '\$E_\\mathrm{ENPT2}\$',
		'out_en' u log(1):4 			w linespoints lc 6 lw 2 title '\$E_\\mathrm{ENPT3}\$',
		2.713319038091934 lw 3 lc 2 title '\$E_\\mathrm{CI}\$'
"
