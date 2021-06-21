set macros;
TMARGIN = "set tmargin at screen 0.90; set bmargin at screen 0.55";
BMARGIN = "set tmargin at screen 0.55; set bmargin at screen 0.20";
LMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.55";
RMARGIN = "set lmargin at screen 0.55; set rmargin at screen 0.95";

set style line 1 lw 2           lc rgb '#C60000'
set style line 2 lw 2           lc rgb '#C65A00'
set style line 3 lw 2           lc rgb '#007777'

#set term epslatex;
#set output 'imgs/gp_aux.eps';

set key spacing 1.5;
set format y '%.1f';
set size 1.0,0.75;

set xlabel 'Iteration (a.\,u.)';
set xrange [-2:35];
set xtics 0,10,100
set yrange [0.5:1.9];
set multiplot layout 1,2;
set size 0.5,0.6;
set key at graph 0.75,0.5;

set label 1 at graph 0.4, 0.20 '$\alpha = 0.0$'

@LMARGIN;
plot 'out_ham_0.0' u 1:2 w linespoints ls 1 title '';
set size 0.5,0.6;
set key at graph 0.9,0.5;

set xrange [-10:100];
set xtics 0,30,160

set label 1 at graph 0.33, 0.23 '\footnotesize{$\alpha = 0.1$}'
set label 2 at graph 0.35, 0.43 '\footnotesize{$\alpha = 0.2$}'
set label 3 at graph 0.40, 0.60 '\footnotesize{$\alpha = 0.3$}'
set label 4 at graph 0.60, 0.73 '\footnotesize{$\alpha = 0.4$}'
set label 5 at graph 0.73, 0.88 '\footnotesize{$\alpha = 0.5$}'

set format y '';
set ylabel '';
@RMARGIN;
#plot 'out_orb_mix' u 1:2 w linespoints ls 2 title 'orb. mix.', \
#     'out_ham_mix' u 1:2 w linespoints ls 3 title 'ham. mix.';
plot \
    'out_orb_0.2' u 1:2 w lines ls 2 lw 3 lc rgb '#666666' title '', \
    'out_orb_0.4' u 1:2 w lines ls 2 lw 3 lc rgb '#222222' title '', \
    'out_orb_0.6' u 1:2 w lines ls 2 lw 4 lc rgb '#000000' title '', \
