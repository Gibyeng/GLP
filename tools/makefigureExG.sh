#!/bin/bash

# exC
gnuplot << EOF
    set terminal pdf
    set output "../KDD-Submission/figures/ExG.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,8" offset 3,0
    set xtics scale 0,0.001
    set xtics font "Time-Roman,6" offset 0.5,0.5
    set ytics font "Time-Roman,6" 
    set logscale y 10
    set key top left
    set key font ",6"
    set size ratio 0.20
    set boxwidth 0.8
    set style fill solid
    plot "../experiment/exG.dat" using (\$2):xtic(1) title "TG" with histogram lt 2 fs pattern 7 ,"../experiment/exG.dat" using (\$4):xtic(1) title "GRAPE" with histogram lt 4 fs pattern 5,"../experiment/exG.dat" using (\$3):xtic(1) title "OMP" with histogram lt 3 fs pattern 9
    quit
EOF