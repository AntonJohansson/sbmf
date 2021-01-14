#!/bin/sh
	#set term epslatex;
	#set output 'imgs/econv_energy.eps';
	#set key at graph 0.9,0.80;
gnuplot -p -e "
	set size 0.75,0.75;
	set key box;
	set key spacing 1.5;
	set key outside right;
	set format y '%.1f';

	set term epslatex;
	set output 'imgs/econv_energy.eps';
	set xlabel 'Basis size (a.\,u.)';
	set ylabel 'Energy per particle (a.\,u.)';
	set yrange [-0.2:*];
	plot
		'out_H_1c_1g' 	u log(1):2:xtic(1) w linespoints lw 2 title 'System 1',
		'out_H_1c_-1g' 	u log(1):2		   w linespoints lw 2 title 'System 2',
		'out_HG_1c_1g' 	u log(1):2		   w linespoints lw 2 title 'System 3',
		'out_HG_1c_-1g' u log(1):2		   w linespoints lw 2 title 'System 4',
		'out_H_2c_+g' 	u log(1):2		   w linespoints lw 2 title 'System 5',
		'out_H_2c_-g' 	u log(1):2		   w linespoints lw 2 title 'System 6',
		'out_HG_2c_+g' 	u log(1):2		   w linespoints lw 2 title 'System 7',
		'out_HG_2c_-g' 	u log(1):2		   w linespoints lw 2 title 'System 8',
"
#set output 'imgs/econv_time.eps';
#set key at graph 0.75,0.85;
#set xlabel 'Basis size (unitless)';
#set ylabel 'Avg. iteration time (sec./iteration)';
#set format y '%.0f';
#plot
#	'out-411-4' u log(1):3:xtic(1) w linespoints lw 2 title 'Strong attr.',
#	'out-211-2' u log(1):3 		   w linespoints lw 2 title 'Weak attr.',
#	'out1221' 	u log(1):3 		   w linespoints lw 2 title 'Weak repul.',
#	'out1441' 	u log(1):3 		   w linespoints lw 2 title 'Strong repul.';
