#!/bin/bash

# exC
gnuplot << EOF
    set terminal pdf
    set output "../Tex paper/figures/dblpExC.pdf"
    set ylabel "Speed up" font "Time-Roman,25"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,15" 
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/dblp/exC.dat" using (\$3/\$2):xtic(1) notitle with boxes lt 2 lc black fs pattern 5
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set output "../Tex paper/figures/livejExC.pdf"
    set ylabel "Speed up" font "Time-Roman,25"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,15" 
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/livej/exC.dat" using (\$3/\$2):xtic(1) notitle with boxes lt 2 lc black fs pattern 5
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set output "../Tex paper/figures/twitterExC.pdf"
    set ylabel "Speed up" font "Time-Roman,25"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,15" 
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/twitter/exC.dat" using (\$3/\$2):xtic(1) notitle with boxes lt 2 lc black fs pattern 5
    quit
EOF

# gnuplot << EOF
#     set terminal pdf
#     set output "../Tex paper/figures/ukExC.pdf"
#     set ylabel "Speed up" font "Time-Roman,25"
#     set xtics font "Time-Roman,25"
#     set ytics font "Time-Roman,15" 
#     set boxwidth 0.3
#     set style fill solid
#     plot "../experiment/uk-2002/exC.dat" using (\$3/\$2):xtic(1) notitle with boxes lt 2 lc black fs pattern 5
#     quit
# EOF


# gnuplot << EOF
#     set terminal pdf
#     set output "../Tex paper/figures/wikiExC.pdf"
#     set ylabel "Speed up" font "Time-Roman,25"
#     set xtics font "Time-Roman,25"
#     set ytics font "Time-Roman,15" 
#     set boxwidth 0.3
#     set style fill solid
#     plot "../experiment/wiki-en/exC.dat" using (\$3/\$2):xtic(1) notitle with boxes lt 2 lc black fs pattern 5
#     quit
# EOF

# gnuplot << EOF
#     set terminal pdf
#     set output "../Tex paper/figures/youtubeExC.pdf"
#     set ylabel "Speed up" font "Time-Roman,25"
#     set xtics font "Time-Roman,25"
#     set ytics font "Time-Roman,15" 
#     set boxwidth 0.3
#     set style fill solid
#     plot "../experiment/youtube/exC.dat" using (\$3/\$2):xtic(1) notitle with boxes lt 2 lc black fs pattern 5
#     quit
# EOF

# gnuplot << EOF
#     set terminal pdf
#     set output "../Tex paper/figures/aligraphExC.pdf"
#     set ylabel "Speed up" font "Time-Roman,25"
#     set xtics font "Time-Roman,25"
#     set ytics font "Time-Roman,15" 
#     set boxwidth 0.3
#     set style fill solid
#     plot "../experiment/aligraph/exC.dat" using (\$3/\$2):xtic(1) notitle with boxes lt 2 lc black fs pattern 5
#     quit
# EOF

# gnuplot << EOF
#     set terminal pdf
#     set output "../Tex paper/figures/roadnetExC.pdf"
#     set ylabel "Speed up" font "Time-Roman,25"
#     set xtics font "Time-Roman,25"
#     set ytics font "Time-Roman,15" 
#     set boxwidth 0.3
#     set style fill solid
#     plot "../experiment/roadnet/exC.dat" using (\$3/\$2):xtic(1) notitle with boxes lt 2 lc black fs pattern 5
#     quit
# EOF