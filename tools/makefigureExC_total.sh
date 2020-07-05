#!/bin/bash
#based on OMP method
# exC
gnuplot << EOF
    set terminal pdf
    set border lw 0.6
    set output "../Sigmod Submission/figures/ExC.pdf"
    set ylabel "Speedup" font "Time-Roman,11" offset 2,0
    set xtics scale 0,0.001
    set xtics font "Time-Roman,11" offset 0.5,0.5
    set ytics font "Time-Roman,11" offset -1,0
    set logscale y
    set yrange [:2000]
    set key top left ins horizontal spacing 0.1
    set key font ",9.1"
    set size ratio 0.20
    set boxwidth 1
    set style fill solid

	plot  "../experiment/exC.dat" using (\$5/\$7):xtic(1) title "Tigergraph" with histogram lt 3 fs pattern 8, "../experiment/exC.dat" using (\$5/\$6):xtic(1) title "Ligra" with histogram lt 6 fs pattern 7, "../experiment/exC.dat" using (\$5/\$2):xtic(1) title "G-Hash" with histogram lt 2 fs pattern 6, "../experiment/exC.dat" using (\$5/\$3):xtic(1) title "G-Sort" with histogram lt 4 fs pattern 9, "../experiment/exC.dat" using (\$5/\$4):xtic(1) title "GLP" with histogram lt 3 fs pattern 3
    quit
EOF


#    set label font "Time-Roman,5" "5.01" at 0,1 offset -2.0,-0.9
#     set label font "Time-Roman,5" "0.66" at 0,1 offset -1.2,0.1
#     set label font "Time-Roman,5" "0.03" at 0,1 offset -0.4,1.8
#     set label font "Time-Roman,5" "0.01" at 0,1 offset 0.6,2.3
#     set label font "Time-Roman,5" "0.003" at 0,1 offset 1.2,3

#     set label font "Time-Roman,5" "27.51" at 1,1 offset -2.4,-1.6
#     set label font "Time-Roman,5" "1.46" at 1,1 offset -1.3,-0.1
#     set label font "Time-Roman,5" "0.15" at 1,1 offset -0.3,1.1
#     set label font "Time-Roman,5" "0.05" at 1,1 offset 0.6,1.7
#     set label font "Time-Roman,5" "0.008" at 1,1 offset 1.4,2.6
    
#     set label font "Time-Roman,5" "16.72" at 2,1 offset -2.4,-1.0
#     set label font "Time-Roman,5" "1.26" at 2,1 offset -1.3,0.4
#     set label font "Time-Roman,5" "0.07" at 2,1 offset -0.3,1.9
#     set label font "Time-Roman,5" "0.03" at 2,1 offset 0.6,2.4
#     set label font "Time-Roman,5" "0.006" at 2,1 offset 1.2,3.1

#     set label font "Time-Roman,5" "19.32" at 3,1 offset -2.4,-0.8
#     set label font "Time-Roman,5" "2.71" at 3,1 offset -1.2,0.2
#     set label font "Time-Roman,5" "0.47" at 3,1 offset -0.3,1.1
#     set label font "Time-Roman,5" "0.24" at 3,1 offset 0.6,1.5
#     set label font "Time-Roman,5" "0.07" at 3,1 offset 1.5,2.2

#     set label font "Time-Roman,5" "157" at 4,1 offset -2.1,-1.0
#     set label font "Time-Roman,5" "18.30" at 4,1 offset -1.6,0.1
#     set label font "Time-Roman,5" "0.37" at 4,1 offset -0.4,2.2
#     set label font "Time-Roman,5" "0.28" at 4,1 offset 0.7,2.5
#     set label font "Time-Roman,5" "0.11" at 4,1 offset 1.6,2.9

#     set label font "Time-Roman,5" "412" at 5,1 offset -2.0,-1.0
#     set label font "Time-Roman,5" "38.52" at 5,1 offset -1.6,0.2
#     set label font "Time-Roman,5" "1.11" at 5,1 offset -0.3,2.0
#     set label font "Time-Roman,5" "0.66" at 5,1 offset 0.7,2.3
#     set label font "Time-Roman,5" "0.23" at 5,1 offset 1.6,2.9

#     set label font "Time-Roman,5" "1343" at 6,1 offset -2.4,-1.1
#     set label font "Time-Roman,5" "102" at 6,1 offset -1.3,0.2
#     set label font "Time-Roman,5" "1.84" at 6,1 offset -0.6,2.4
#     set label font "Time-Roman,5" "1.68" at 6,1 offset 0.6,2.5
#     set label font "Time-Roman,5" "0.84" at 6,1 offset 1.5,2.9

#     set label font "Time-Roman,5" "3409" at 7,1 offset -2.4,-1.3
#     set label font "Time-Roman,5" "202" at 7,1 offset -1.2,0.2
#     set label font "Time-Roman,5" "8.30" at 7,1 offset -0.5,1.9
#     set label font "Time-Roman,5" "10.1" at 7,1 offset 0.7,1.8
#     set label font "Time-Roman,5" "1.78" at 7,1 offset 1.5,2.7