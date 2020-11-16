#!/bin/sh
gnuplot -p -e "
	plot
	'out1' u 1:2,
	'out2' u 1:2,
	'out3' u 1:2,
	'out4' u 1:2,
	'out5' u 1:2,
	'out6' u 1:2
"
