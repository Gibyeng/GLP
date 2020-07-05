#!/bin/bash
gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/dblpExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "pull push fg-push ts-push"
    set xlabel font ",18" offset 2.5, 0
    plot "../experiment/dblp/exAplus.dat" using 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, '' u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 12, ''  u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 5, '' u 4 with histograms lt 2 lw 2 lc rgb "red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/livejExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "pull push fg-push ts-push"
    set xlabel font ",18" offset 2.1, 0
    plot "../experiment/livej/exAplus.dat" using 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, '' u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 12, ''  u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 5, '' u 4 with histograms lt 2 lw 2 lc rgb "red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/twitterExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "pull push fg-push ts-push"
    set xlabel font ",18" offset 2.1, 0
    plot "../experiment/twitter/exAplus.dat" using 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, '' u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 12, ''  u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 5, '' u 4 with histograms lt 2 lw 2 lc rgb "red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/ukExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "pull push fg-push ts-push"
    set xlabel font ",18" offset 2.1, 0
    plot "../experiment/uk-2002/exAplus.dat" using 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, '' u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 12, ''  u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 5, '' u 4 with histograms lt 2 lw 2 lc rgb "red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/wikiExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "pull push fg-push ts-push"
    set xlabel font ",18" offset 2.1, 0
    plot "../experiment/wiki-en/exAplus.dat" using 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, '' u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 12, ''  u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 5, '' u 4 with histograms lt 2 lw 2 lc rgb "red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/youtubeExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "pull push fg-push ts-push"
    set xlabel font ",18" offset 2.1, 0
    plot "../experiment/youtube/exAplus.dat" using 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, '' u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 12, ''  u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 5, '' u 4 with histograms lt 2 lw 2 lc rgb "red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/aligraphExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "pull push fg-push ts-push"
    set xlabel font ",18" offset 2.1, 0
    plot "../experiment/aligraph/exAplus.dat" using 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, '' u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 12, ''  u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 5, '' u 4 with histograms lt 2 lw 2 lc rgb "red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/roadnetExA.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "pull push fg-push ts-push"
    set xlabel font ",18" offset 2.1, 0
    plot "../experiment/roadnet/exAplus.dat" using 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, '' u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 12, ''  u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 5, '' u 4 with histograms lt 2 lw 2 lc rgb "red" fs pattern 1
    quit
EOF