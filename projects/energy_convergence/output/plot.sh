#!/bin/sh
gnuplot -p -e "
	set term epslatex;
	set output 'imgs/econv_energy.eps';
	plot
		'out-411-4' u 1:2:xtic(1) w linespoints title 'strong attr.',
		'out-211-2' u 1:2 		  w linespoints title 'weak attr.',
		'out1221' 	u 1:2 		  w linespoints title 'weak repul.',
		'out1441' 	u 1:2 		  w linespoints title 'strong repul.',
"
gnuplot -p -e "
	set term epslatex;
	set output 'imgs/econv_time.eps';
	plot
		'out-411-4' u 1:3:xtic(1) w linespoints title 'strong attr.',
		'out-211-2' u 1:3 		  w linespoints title 'weak attr.',
		'out1221' 	u 1:3 		  w linespoints title 'weak repul.',
		'out1441' 	u 1:3 		  w linespoints title 'strong repul.',
"
