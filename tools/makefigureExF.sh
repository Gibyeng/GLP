#!/bin/bash

# exF
gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/dblpExF.pdf"
    set ylabel "Reduced memory access(%)" font "Time-Roman,18"
    set xlabel "d of CMS" font "Time-Roman,21"
    set yrange [84:96]
    set xtics font "Time-Roman,18"
    set ytics font "Time-Roman,18" 
    set key font ",14"
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/dblp/exF.dat" using (100-(\$2/\$5*100)):xtic(1) title "r=0.5" with lp lw 3 lt 2 ps 1.5, "../experiment/dblp/exF.dat" using (100-(\$3/\$5*100)) title "r=   1" with lp lc 'red' lw 3 lt 8 ps 1.5, "../experiment/dblp/exF.dat" using (100-(\$4/\$5*100)) title "r=   2" with lp lw 3 lt 4 ps 1.5   
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/livejExF.pdf"
    set ylabel "Reduced memory access(%)" font "Time-Roman,18"
    set xlabel "d of CMS" font "Time-Roman,21"
    set yrange [:86]
    set xtics font "Time-Roman,18"
    set ytics font "Time-Roman,18" 
    set key font ",14"
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/livej/exF.dat" using (100-(\$2/\$5*100)):xtic(1) title "r=0.5" with lp lw 3 lt 2 ps 1.5, "../experiment/livej/exF.dat" using (100-(\$3/\$5*100)) title "r=   1" with lp lc 'red' lw 3 lt 8 ps 1.5, "../experiment/livej/exF.dat" using (100-(\$4/\$5*100)) title "r=   2" with lp lw 3 lt 4 ps 1.5   
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/twitterExF.pdf"
    set ylabel "Reduced memory access(%)" font "Time-Roman,18"
    set xlabel "d of CMS" font "Time-Roman,21"
    set yrange [80:90]
    set xtics font "Time-Roman,18"
    set ytics font "Time-Roman,18"
    set key font ",14" 
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/twitter/exF.dat" using (100-(\$2/\$5*100)):xtic(1) title "r=0.5" with lp lw 3 lt 2 ps 1.5, "../experiment/twitter/exF.dat" using (100-(\$3/\$5*100)) title "r=   1" with lp lc 'red' lw 3 lt 8 ps 1.5, "../experiment/twitter/exF.dat" using (100-(\$4/\$5*100)) title "r=   2" with lp lw 3 lt 4 ps 1.5   
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/ukExF.pdf"
    set ylabel "Reduced memory access(%)" font "Time-Roman,18"
    set xlabel "d of CMS" font "Time-Roman,21"
    set yrange [84:94]
    set xtics font "Time-Roman,18"
    set ytics font "Time-Roman,18"
    set key font ",14" 
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/uk-2002/exF.dat" using (100-(\$2/\$5*100)):xtic(1) title "r=0.5" with lp lw 3 lt 2 ps 1.5, "../experiment/uk-2002/exF.dat" using (100-(\$3/\$5*100)) title "r=   1" with lp lc 'red' lw 3 lt 8 ps 1.5, "../experiment/uk-2002/exF.dat" using (100-(\$4/\$5*100)) title "r=   2" with lp lw 3 lt 4 ps 1.5   
    quit
EOF


gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/wikiExF.pdf"
    set ylabel "Reduced memory access(%)" font "Time-Roman,18"
    set xlabel "d of CMS" font "Time-Roman,21"
    set yrange [:*]
    set xtics font "Time-Roman,18"
    set ytics font "Time-Roman,18"
    set key font ",14" 
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/wiki-en/exF.dat" using (100-(\$2/\$5*100)):xtic(1) title "r=0.5" with lp lw 3 lt 2 ps 1.5, "../experiment/wiki-en/exF.dat" using (100-(\$3/\$5*100)) title "r=   1" with lp lc 'red' lw 3 lt 8 ps 1.5, "../experiment/wiki-en/exF.dat" using (100-(\$4/\$5*100)) title "r=   2" with lp lw 3 lt 4 ps 1.5    
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/youtubeExF.pdf"
    set ylabel "Reduced memory access(%)" font "Time-Roman,18"
    set xlabel "d of CMS" font "Time-Roman,21"
    set yrange [70:86]
    set xtics font "Time-Roman,18"
    set ytics font "Time-Roman,18"
    set key font ",14" 
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/youtube/exF.dat" using (100-(\$2/\$5*100)):xtic(1) title "r=0.5" with lp lw 3 lt 2 ps 1.5, "../experiment/youtube/exF.dat" using (100-(\$3/\$5*100)) title "r=   1" with lp lc 'red' lw 3 lt 8 ps 1.5, "../experiment/youtube/exF.dat" using (100-(\$4/\$5*100)) title "r=   2" with lp lw 3 lt 4 ps 1.5    
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/roadnetExF.pdf"
    set ylabel "Reduced memory access(%)" font "Time-Roman,18"
    set xlabel "d of CMS" font "Time-Roman,21"
    set yrange [:*]
    set xtics font "Time-Roman,18"
    set ytics font "Time-Roman,18"
    set key font ",14" 
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/roadnet/exF.dat" using (100-(\$2/\$5*100)):xtic(1) title "r=0.5" with lp lw 3 lt 2 ps 1.5, "../experiment/roadnet/exF.dat" using (100-(\$3/\$5*100)) title "r=   1" with lp lc 'red' lw 3 lt 8 ps 1.5, "../experiment/roadnet/exF.dat" using (100-(\$4/\$5*100)) title "r=   2" with lp lw 3 lt 4 ps 1.5    
    quit
EOF

gnuplot << EOF
    set terminal pdf
    set border lw 3
    set output "../KDD-Submission/figures/aligraphExF.pdf"
    set ylabel "Reduced memory access(%)" font "Time-Roman,18"
    set xlabel "d of CMS" font "Time-Roman,21"
    set yrange [:*]
    set xtics font "Time-Roman,18"
    set ytics font "Time-Roman,18"
    set key font ",14" 
    set boxwidth 0.3
    set style fill solid
    plot "../experiment/aligraph/exF.dat" using (100-(\$2/\$5*100)):xtic(1) title "r=0.5" with lp lw 3 lt 2 ps 1.5, "../experiment/aligraph/exF.dat" using (100-(\$3/\$5*100)) title "r=   1" with lp lc 'red' lw 3 lt 8 ps 1.5, "../experiment/aligraph/exF.dat" using (100-(\$4/\$5*100)) title "r=   2" with lp lw 3 lt 4 ps 1.5  
    quit
EOF
