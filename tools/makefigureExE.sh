#!/bin/bash



gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/ExE_m_16.pdf"
    set xlabel "Number of vertices(log)" font "Time-Roman,18"
    set ylabel "Elapsed time(s)" font "Time-Roman,25"
    set xtics font "Time-Roman,18"
    set ytics font "Time-Roman,18"
    set xtics nomirror
    set logscale y 
    set boxwidth 0.5
    set key left top ins horizontal spacing 0.7
    set key width 5
    set key font ",13"
    set style fill solid
    plot "../experiment/scalability/exE_m_16.dat" using 3:xtic(1) title 'OMP' with histograms lt 6 lw 2 fs pattern 5, "../experiment/scalability/exE_m_16.dat" using 4 title 'G-Hash' with histograms lt 2 lw 2 fs pattern 6, "../experiment/scalability/exE_m_16.dat" using 5 title 'GLP' with histograms lt 3 lw 2 fs pattern 3, "../experiment/scalability/exE_m_16.dat" using 6 title 'GPU2GPU' with histograms lt 7 lw 2 fs pattern 8
    quit
EOF



gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/ExE_m.pdf"
    set xlabel "Number of edges(hundred million)x" font "Time-Roman,18"
    set ylabel "Elapsed time(s)" font "Time-Roman,25"
    set xtics font "Time-Roman,18"
    set ytics font "Time-Roman,18" 
    set xtics nomirror
    set logscale y
    set boxwidth 0.5
    set key left top ins horizontal spacing 0.7
    set key width 5
    set key font ",13"
    set style fill solid
    plot "../experiment/scalability/exE_m.dat" using 3:xtic(1) title 'OMP' with histograms lt 6 lw 2 fs pattern 5,"../experiment/scalability/exE_m.dat" using 4:xtic(1) title 'G-Hash' with histograms lt 2 lw 2 fs pattern 6 , "../experiment/scalability/exE_m.dat" using 5 title 'GLP' with histograms lt 3 lw 2 fs pattern 3,"../experiment/scalability/exE_m.dat" using 6 title 'GPU2GPU' with histograms lt 7 lw 2 fs pattern 8
    quit
EOF

