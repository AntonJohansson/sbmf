#!/bin/sh

	#set term epslatex;
	#set output 'imgs/gp_pt_particle_conv.eps';
gnuplot -p -e "
	set size 0.75,0.75;
	set key at graph 0.9,0.90;
	set key box;
	set format y '%.3e';
	set grid;

	set xlabel 'Basis size (a.\,u.)';
	set ylabel 'Energy (a.\,u.)';
	plot
		'out_32' u 1:4 	w linespoints lw 2 title '2',
"
