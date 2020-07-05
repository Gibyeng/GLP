#!/bin/bash

# exOpt
gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../Sigmod Submission/figures/dblpExOpt.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/dblp/exOpt.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../Sigmod Submission/figures/livejExOpt.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/livej/exOpt.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../Sigmod Submission/figures/twitterExOpt.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" 
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/twitter/exOpt.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../Sigmod Submission/figures/ukExOpt.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" 
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/uk-2002/exOpt.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF


gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../Sigmod Submission/figures/wikiExOpt.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" 
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/wiki-en/exOpt.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../Sigmod Submission/figures/youtubeExOpt.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/youtube/exOpt.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../Sigmod Submission/figures/aligraphExOpt.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/aligraph/exOpt.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../Sigmod Submission/figures/roadnetExOpt.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,20"
    set ytics font "Time-Roman,15" 
    set yrange [0:*]
    set boxwidth 0.2
    set style fill solid
    plot "../experiment/roadnet/exOpt.dat" using (\$0):2:3:xtic(1) notitle with boxes lt 2 lw 2 lc variable fs pattern 3
    quit
EOF