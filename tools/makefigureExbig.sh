#!/bin/bash

# exC
gnuplot << EOF
    set terminal pdf
    set border lw 0.8
    set output "../Sigmod Submission/figures/Exbig.pdf"
    set xlabel "Days" font "Time-Roman,11" offset 0,1
    set ylabel "Elapsed time(s)" font "Time-Roman,11" offset 2,0
    set xtics scale 0,0.001
    set xtics font "Time-Roman,11" offset 0.5,0.5
    set ytics font "Time-Roman,11" offset 0,0
    set logscale y
    set key top left ins horizontal spacing 0.1
    set key font ",9.1"
    set size ratio 0.20
    set boxwidth 0.8
    set style fill solid
    
	plot "../experiment/exbig.dat" using (\$2):xtic(1) title "GLP" with histogram lt 3 fs pattern 3,"../experiment/exbig.dat" using (\$4):xtic(1) title "ALI" with histogram lt 10 fs pattern 4
    quit
EOF

    # set label font "Time-Roman,8" "10.0%" at 0,1 offset -0.8,0.2
    # set label font "Time-Roman,8" "9.8%" at 1,1 offset -0.8,0.6
    # set label font "Time-Roman,8" "9.8%" at 2,1 offset -0.8,0.9
    # set label font "Time-Roman,8" "9.8%" at 3,1 offset -0.8,1.1
    # set label font "Time-Roman,8" "9.9%" at 4,1 offset -0.8,1.3
    # set label font "Time-Roman,8" "9.9%" at 5,1 offset -0.8,1.4
    # set label font "Time-Roman,8" "10.0%" at 6,1 offset -0.8,1.5
    # set label font "Time-Roman,8" "10.0%" at 7,1 offset -0.8,1.6
    # set label font "Time-Roman,8" "10.0%" at 8,1 offset -0.8,1.7
    # set label font "Time-Roman,8" "10.0%" at 9,1 offset -0.8,1.9