#!/bin/bash

# exA
gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/dblpExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/dblp/exA.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/livejExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/livej/exA.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/twitterExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/twitter/exA.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/ukExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15"
    set yrange [0:*] 
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/uk-2002/exA.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF


gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/wikiExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/wiki-en/exA.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/youtubeExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,24"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/youtube/exA.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/aligraphExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/aligraph/exA.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/roadnetExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,24"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/roadnet/exA.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF
