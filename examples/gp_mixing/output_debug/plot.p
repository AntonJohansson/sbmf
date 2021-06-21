set macros;
TMARGIN = "set tmargin at screen 0.90; set bmargin at screen 0.55";
BMARGIN = "set tmargin at screen 0.55; set bmargin at screen 0.20";
TBMARGIN = "set tmargin at screen 0.70; set bmargin at screen 0.20";
#LMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.55";
#RMARGIN = "set lmargin at screen 0.55; set rmargin at screen 0.95";

LMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.4166";
MMARGIN = "set lmargin at screen 0.4166; set rmargin at screen 0.6833";
RMARGIN = "set lmargin at screen 0.6833; set rmargin at screen 0.95";

set style line 1 lw 4           lc rgb '#582A15'
set style line 2 lw 3           lc rgb '#C65A00'
set style line 3 lw 3           lc rgb '#007777'

set term epslatex;
set output 'imgs/gp_aux.eps';
set key spacing 1.5;
set format y '%.1f';
set size 1.0,0.75;

set xlabel 'Iteration';
set ylabel 'Eigenvalue $\mu$ (a.\,u.)';
set yrange [1.2:1.9];

set multiplot layout 1,3;

set label 1 at graph 0.25, 0.20 'pure SCF'

set xrange [-2:25];
set xtics 0,10,30

@TBMARGIN; @LMARGIN;
plot 'out_ham_0.0' u 1:2 w linespoints ls 1 title '';
set key at graph 0.9,0.5;

set xrange [-10:170];
set xtics 0,75,150

set label 1 at graph 0.33, 0.23 '\footnotesize{$\alpha = 0.1$}'
set label 2 at graph 0.35, 0.43 '\footnotesize{$\alpha = 0.2$}'
set label 3 at graph 0.40, 0.60 '\footnotesize{$\alpha = 0.3$}'
set label 4 at graph 0.60, 0.73 '\footnotesize{$\alpha = 0.4$}'
set label 5 at graph 0.60, 0.88 '\footnotesize{$\alpha = 0.5$}'

set format y '';
set ylabel '';
@TBMARGIN; @MMARGIN;
plot \
    'out_ham_0.1' u 1:2 w lines ls 2 lw 3 lc rgb '#FFA19C' title '', \
    'out_ham_0.2' u 1:2 w lines ls 2 lw 3 lc rgb '#EE706A' title '', \
    'out_ham_0.3' u 1:2 w lines ls 2 lw 3 lc rgb '#CE4841' title '', \
    'out_ham_0.4' u 1:2 w lines ls 2 lw 3 lc rgb '#AF2922' title '', \
    'out_ham_0.5' u 1:2 w lines ls 2 lw 4 lc rgb '#890F08' title '', \

#set label 1 front at graph 0.115, 0.10 '\scriptsize{$\beta = 0.90$}'
#set label 2 front at graph 0.455, 0.60 '\scriptsize{$\beta = 0.95$}'
unset label 1
unset label 2
unset label 3
unset label 4
unset label 5

#set arrow 1 front from graph 0.65,0.1 to graph 0.78,0.1 lw 2 lc rgb '#4E9D5C'
#set arrow 2 front from graph 0.65,0.2 to graph 0.82,0.3 lw 2 lc rgb '#01400C'

set xrange [-100:950]
set xtics 0,400,1000
set key spacing 1
set key width -1
set key at graph 0.95,0.25 box opaque samplen 2
set key invert
@TBMARGIN; @RMARGIN;
plot \
    'out_orb_0.9'  u 1:2 w lines ls 2 lw 3 lc rgb '#4E9D5C' title '\scriptsize{$\beta=0.90$}', \
    'out_orb_0.95' u 1:2 w lines ls 2 lw 4 lc rgb '#01400C' title '\scriptsize{$\beta=0.95$}', \
