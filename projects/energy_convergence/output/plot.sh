#!/bin/sh
	#set term epslatex;
	#set output 'imgs/econv_energy.eps';
gnuplot -p -e "
	set term epslatex;
	set output 'imgs/econv_energy.eps';
	set size 0.75,0.75;
	set key at graph 0.9,0.80;
	set key box;
	set grid;

	set output 'imgs/econv_energy.eps';
	set xlabel 'Basis size (unitless)';
	set ylabel 'Energy per particle (a.u.)';
	set yrange [-0.2:*];
	plot
		'out-411-4' u log(1):2:xtic(1) w linespoints lw 2 title 'strong attr.',
		'out-211-2' u log(1):2 		   w linespoints lw 2 title 'weak attr.',
		'out1221' 	u log(1):2 		   w linespoints lw 2 title 'weak repul.',
		'out1441' 	u log(1):2 		   w linespoints lw 2 title 'strong repul.';

	set output 'imgs/econv_time.eps';
	set key at graph 0.75,0.85;
	set xlabel 'Basis size (unitless)';
	set ylabel 'Avg. iteration time (sec./iteration)';
	plot
		'out-411-4' u log(1):3:xtic(1) w linespoints lw 2 title 'strong attr.',
		'out-211-2' u log(1):3 		   w linespoints lw 2 title 'weak attr.',
		'out1221' 	u log(1):3 		   w linespoints lw 2 title 'weak repul.',
		'out1441' 	u log(1):3 		   w linespoints lw 2 title 'strong repul.';

"
