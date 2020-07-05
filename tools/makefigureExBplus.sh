#!/bin/bash
gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/dblpExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "global  smem warp-centric"
    set xlabel font ",18" offset 3, 0
    plot "../experiment/dblp/exBplus.dat" u 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, ''  u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 5, '' u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/livejExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "global  smem warp-centric"
    set xlabel font ",18" offset 3, 0
    plot "../experiment/livej/exBplus.dat" u 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, ''  u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 5, '' u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/twitterExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "global  smem warp-centric"
    set xlabel font ",18" offset 3, 0
    plot "../experiment/twitter/exBplus.dat" u 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, ''  u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 5, '' u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/ukExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "global  smem warp-centric"
    set xlabel font ",18" offset 3, 0
    plot "../experiment/uk-2002/exBplus.dat" u 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, ''  u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 5, '' u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/wikiExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "global  smem warp-centric"
    set xlabel font ",18" offset 3, 0
    plot "../experiment/wiki-en/exBplus.dat" u 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, ''  u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 5, '' u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/youtubeExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "global  smem warp-centric"
    set xlabel font ",18" offset 3, 0
    plot "../experiment/youtube/exBplus.dat" u 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, ''  u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 5, '' u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/aligraphExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "global  smem warp-centric"
    set xlabel font ",18" offset 3, 0
    plot "../experiment/aligraph/exBplus.dat" u 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, ''  u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 5, '' u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 1
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set size 1,1
    set output "../KDD-Submission/figures/roadnetExB.pdf"
    set ylabel "Elapsed time(s)" font "Time-Roman,23" offset "-0.6,0"
    set xtics font "Time-Roman,25"
    set ytics font "Time-Roman,18" 
    set yrange [0:*]
    set boxwidth 0.6
    set style histogram clustered gap 0
    unset xtics
    unset key
    set xlabel "global  smem warp-centric"
    set xlabel font ",18" offset 3, 0
    plot "../experiment/roadnet/exBplus.dat" u 1 with histograms lt 2 lw 2 lc rgb "tan1" fs pattern 8, ''  u 2 with histograms lt 2 lw 2 lc rgb "salmon" fs pattern 5, '' u 3 with histograms lt 2 lw 2 lc rgb "orange-red" fs pattern 1
    quit
EOF