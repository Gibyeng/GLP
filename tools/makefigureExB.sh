#!/bin/bash

# exB
gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/dblpExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/dblp/exB.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/livejExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/livej/exB.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/twitterExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" 
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/twitter/exB.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/ukExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" 
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/uk-2002/exB.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF


gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/wikiExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" 
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/wiki-en/exB.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/youtubeExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/youtube/exB.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/aligraphExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/aligraph/exB.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/roadnetExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/roadnet/exB.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF