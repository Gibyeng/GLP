#!/bin/bash
#based on OMP method
# exlayered
gnuplot << EOF
    set terminal pdf
    set border lw 0.6
    set output "../Sigmod Submission/figures/Exlayered.pdf"
    set ylabel "Speedup" font "Time-Roman,11" offset 2,0
    set xtics scale 0,0.001
    set xtics font "Time-Roman,11" offset 0.5,0.5
    set ytics font "Time-Roman,11" offset -1,0
    set logscale y
    set key top left ins horizontal spacing 0.1
    set key font ",9.1"
    set size ratio 0.20
    set boxwidth 1
    set style fill solid

	plot "../experiment/exlayered.dat" using (\$5/\$6):xtic(1) title "Ligra" with histogram lt 6 fs pattern 7, "../experiment/exlayered.dat" using (\$5/\$2):xtic(1) title "G-Hash" with histogram lt 2 fs pattern 6, "../experiment/exlayered.dat" using (\$5/\$3):xtic(1) title "G-Sort" with histogram lt 4 fs pattern 9, "../experiment/exlayered.dat" using (\$5/\$4):xtic(1) title "GLP" with histogram lt 3 fs pattern 3
    quit
EOF

#     set label font "Time-Roman,5" "5.97" at 0,1 offset -1.8,0.2
#     set label font "Time-Roman,5" "0.24" at 0,1 offset -0.8,2.4
#     set label font "Time-Roman,5" "0.17" at 0,1 offset 0.4,2.7
#     set label font "Time-Roman,5" "0.05" at 0,1 offset 1.5,3.5

#     set label font "Time-Roman,5" "14.13" at 1,1 offset -2,0.2
#     set label font "Time-Roman,5" "1.10" at 1,1 offset -0.7,2.0
#     set label font "Time-Roman,5" "0.78" at 1,1 offset 0.5,2.2
#     set label font "Time-Roman,5" "0.19" at 1,1 offset 1.4,3.2
    
#     set label font "Time-Roman,5" "22.26" at 2,1 offset -2,0
#     set label font "Time-Roman,5" "0.66" at 2,1 offset -0.8,2.4
#     set label font "Time-Roman,5" "0.38" at 2,1 offset 0.3,2.8
#     set label font "Time-Roman,5" "0.22" at 2,1 offset 1.4,3.2

#     set label font "Time-Roman,5" "23.35" at 3,1 offset -2,0.3
#     set label font "Time-Roman,5" "4.44" at 3,1 offset -0.7,1.4
#     set label font "Time-Roman,5" "2.57" at 3,1 offset 0.4,1.8
#     set label font "Time-Roman,5" "0.67" at 3,1 offset 1.5,2.7

#     set label font "Time-Roman,5" "205" at 4,1 offset -1.6,0.1
#     set label font "Time-Roman,5" "3.36" at 4,1 offset -0.8,3
#     set label font "Time-Roman,5" "3.51" at 4,1 offset 0.5,3
#     set label font "Time-Roman,5" "1.30" at 4,1 offset 1.6,3.6

#     set label font "Time-Roman,5" "442" at 5,1 offset -1.6,0.2
#     set label font "Time-Roman,5" "9.70" at 5,1 offset -0.7,2.8
#     set label font "Time-Roman,5" "8.66" at 5,1 offset 0.5,2.9
#     set label font "Time-Roman,5" "3.45" at 5,1 offset 1.6,3.5

#     set label font "Time-Roman,5" "942" at 6,1 offset -1.6,0.2
#     set label font "Time-Roman,5" "16.94" at 6,1 offset -1.1,3
#     set label font "Time-Roman,5" "19.90" at 6,1 offset 0.4,2.9
#     set label font "Time-Roman,5" "9.05" at 6,1 offset 1.5,3.4

#     set label font "Time-Roman,5" "1901" at 7,1 offset -1.8,0.1
#     set label font "Time-Roman,5" "83.15" at 7,1 offset -0.8,2.3
#     set label font "Time-Roman,5" "117" at 7,1 offset 0.6,2.1
#     set label font "Time-Roman,5" "18.67" at 7,1 offset 1.5,3.4
