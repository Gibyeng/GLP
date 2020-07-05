#!/bin/bash
gnuplot << EOF
    set terminal pdf
    set size 0.8,1
    set output "./1.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "global    smem warp-centric"
    set xlabel font ",11" offset 1.8, 0
    plot "./1.txt" u 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, ''  u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 5, '' u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 1
    quit
EOF