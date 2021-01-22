#!/bin/sh
	#set term epslatex;
	#set output 'imgs/econv_energy.eps';
	#set key at graph 0.9,0.80;
gnuplot -p -e "
	set term epslatex;
	set output 'imgs/econv_energy.eps';
	set size 0.75,0.75;
	set key box;
	set key spacing 1.5;
	set key outside right;
	set format y '%.1f';

	set xlabel 'Basis size (a.\,u.)';
	set ylabel 'Energy per particle (a.\,u.)';
	set yrange [0:1.7];
	plot
		'out_1_1c_VH' 	u log(1):2:xtic(1) w linespoints lw 2 title 'Sys. 1',
		'out_-1_1c_VH' 	u log(1):2		   w linespoints lw 2 title 'Sys. 2',
		'out_1_1c_VHG'	u log(1):2		   w linespoints lw 2 title 'Sys. 3',
		'out_-1_1c_VHG'	u log(1):2		   w linespoints lw 2 title 'Sys. 4',
		'out_+_2c_VH'	u log(1):2		   w linespoints lw 2 title 'Sys. 5',
		'out_-_2c_VH'	u log(1):2		   w linespoints lw 2 title 'Sys. 6',
		'out_+_2c_VHG'	u log(1):2		   w linespoints lw 2 title 'Sys. 7',
		'out_-_2c_VHG'	u log(1):2		   w linespoints lw 2 title 'Sys. 8',
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
