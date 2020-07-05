#!/bin/bash
#based on OMP method
# exslpa
gnuplot << EOF
    set terminal pdf
    set border lw 0.6
    set output "../Sigmod Submission/figures/Exslpa.pdf"
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

	plot "../experiment/exslpa.dat" using (\$5/\$6):xtic(1) title "Ligra" with histogram lt 6 fs pattern 7, "../experiment/exslpa.dat" using (\$5/\$2):xtic(1) title "G-Hash" with histogram lt 2 fs pattern 6, "../experiment/exslpa.dat" using (\$5/\$3):xtic(1) title "G-Sort" with histogram lt 4 fs pattern 9, "../experiment/exslpa.dat" using (\$5/\$4):xtic(1) title "GLP" with histogram lt 3 fs pattern 3
    quit
EOF

    # set label font "Time-Roman,5" "0.87" at 0,1 offset -1.8,0.1
    # set label font "Time-Roman,5" "0.03" at 0,1 offset -0.7,2.4
    # set label font "Time-Roman,5" "0.01" at 0,1 offset 0.5,2.9
    # set label font "Time-Roman,5" "0.004" at 0,1 offset 1.2,3.9

    # set label font "Time-Roman,5" "1.48" at 1,1 offset -1.7,-0.2
    # set label font "Time-Roman,5" "0.15" at 1,1 offset -0.7,1.4
    # set label font "Time-Roman,5" "0.05" at 1,1 offset 0.5,2.1
    # set label font "Time-Roman,5" "0.02" at 1,1 offset 1.5,2.9
    
    # set label font "Time-Roman,5" "1.39" at 2,1 offset -1.7,0.4
    # set label font "Time-Roman,5" "0.09" at 2,1 offset -0.7,2.3
    # set label font "Time-Roman,5" "0.03" at 2,1 offset 0.4,3
    # set label font "Time-Roman,5" "0.01" at 2,1 offset 1.4,3.7

    # set label font "Time-Roman,5" "3.45" at 3,1 offset -1.6,0.2
    # set label font "Time-Roman,5" "0.50" at 3,1 offset -0.7,1.5
    # set label font "Time-Roman,5" "0.28" at 3,1 offset 0.4,1.9
    # set label font "Time-Roman,5" "0.08" at 3,1 offset 1.5,2.9

    # set label font "Time-Roman,5" "24.52" at 4,1 offset -1.9,-0.1
    # set label font "Time-Roman,5" "0.43" at 4,1 offset -0.8,2.8
    # set label font "Time-Roman,5" "0.37" at 4,1 offset 0.4,3
    # set label font "Time-Roman,5" "0.13" at 4,1 offset 1.6,3.6

    # set label font "Time-Roman,5" "41.24" at 5,1 offset -2,0.2
    # set label font "Time-Roman,5" "1.63" at 5,1 offset -0.6,2.4
    # set label font "Time-Roman,5" "0.97" at 5,1 offset 0.5,2.8
    # set label font "Time-Roman,5" "0.28" at 5,1 offset 1.6,3.7

    # set label font "Time-Roman,5" "110" at 6,1 offset -1.7,0.2
    # set label font "Time-Roman,5" "2.56" at 6,1 offset -0.7,2.8
    # set label font "Time-Roman,5" "2.10" at 6,1 offset 0.5,3
    # set label font "Time-Roman,5" "0.87" at 6,1 offset 1.5,3.5

    # set label font "Time-Roman,5" "218" at 7,1 offset -1.6,0.1
    # set label font "Time-Roman,5" "9.77" at 7,1 offset -0.9,2.4
    # set label font "Time-Roman,5" "10.61" at 7,1 offset 0.3,2.3
    # set label font "Time-Roman,5" "1.82" at 7,1 offset 1.5,3.4