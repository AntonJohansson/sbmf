#!/bin/sh

gnuplot -p -e "set title 'Memory Usage'; plot 'memory.log' u 1:2 w lines, 'memory.log' u 1:3 w lines"
