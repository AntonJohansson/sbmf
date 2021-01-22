#!/bin/sh
gnuplot -p -e "
	set term epslatex;
	set output 'imgs/gp_aux.eps';
	set key spacing 1.5;
	set format y '%.1f';
	set size 1.0,0.6;

	set xlabel 'Iteration (a.\,u.)';
	set ylabel 'Chemical potential \$\\mu\$ (a.\,u.)';
	set xrange [0:30];
	set yrange [0.5:2];
	set multiplot layout 1,2;
	set size 0.5,0.6;
	set key at graph 0.75,0.5;
	plot
		'out_no_aux' 	u 1:2 w linespoints lw 2 title 'SCF';
	set size 0.5,0.6;
	set key at graph 0.9,0.5;
	plot
		'out_diis_2' 	u 1:2 w linespoints lw 2 lc 2 title 'DIIS',
		'out_orb_mix' 	u 1:2 w linespoints lw 2 lc 4 title 'orb. mix.',
		'out_ham_mix' 	u 1:2 w linespoints lw 2 lc 7 title 'ham. mix.';
"
