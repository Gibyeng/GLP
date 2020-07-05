#!/bin/bash

# exD
gnuplot << EOF
    set terminal pdf
    set output "../KDD-Submission/figures/dblpExD.pdf"
    set ylabel "Accessed Global Memory(byte)" font "Time-Roman,20"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/dblp/exD.dat" using (\$2*\$3):xtic(1) notitle with boxes lt 2
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set output "../KDD-Submission/figures/livejExD.pdf"
    set ylabel "Accessed Global Memory(byte)" font "Time-Roman,20"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/livej/exD.dat" using (\$2*\$3):xtic(1) notitle with boxes lt 2
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set output "../KDD-Submission/figures/twitterExD.pdf"
    set ylabel "Accessed Global Memory(byte)" font "Time-Roman,20"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/twitter/exD.dat" using (\$2*\$3):xtic(1) notitle with boxes lt 2
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set output "../KDD-Submission/figures/ukExD.pdf"
    set ylabel "Accessed Global Memory(byte)" font "Time-Roman,20"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/uk-2002/exD.dat" using (\$2*\$3):xtic(1) notitle with boxes lt 2
    quit
EOF


gnuplot << EOF
    set terminal pdf
    set output "../KDD-Submission/figures/wikiExD.pdf"
    set ylabel "Accessed Global Memory(byte)" font "Time-Roman,20"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/wiki-en/exD.dat" using (\$2*\$3):xtic(1) notitle with boxes lt 2
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set output "../KDD-Submission/figures/youtubeExD.pdf"
    set ylabel "Accessed Global Memory(byte)" font "Time-Roman,20"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/youtube/exD.dat" using (\$2*\$3):xtic(1) notitle with boxes lt 2
    quit
EOF
