#!/bin/bash
# -*- coding: utf-8 -*-

./graphGenerator D1_part1.txt 100000000 400000000 1 1
./graphGenerator D1_part2.txt 100000000 400000000 1 2
./graphGenerator D1_part3.txt 100000000 400000000 1 3
./graphGenerator D1_part4.txt 160000000 500000000 1 4
./graphGenerator D2_part1.txt 50000000 400000000 1 5
./graphGenerator D2_part2.txt 50000000 400000000 1 6
./graphGenerator D2_part3.txt 70000000 500000000 1 7
./graphGenerator D3_part1.txt 20000000 400000000 1 8
./graphGenerator D3_part2.txt 20000000 400000000 1 9
./graphGenerator D3_part3.txt 30000000 500000000 1 10
./graphGenerator D4_part1.txt 20000000 400000000 1 11
./graphGenerator D4_part2.txt 20000000 400000000 1 12
./graphGenerator D4_part3.txt 30000000 400000000 1 13
./graphGenerator D5_part1.txt 20000000 400000000 1 14
./graphGenerator D5_part2.txt 20000000 400000000 1 15
./graphGenerator D5_part3.txt 10000000 400000000 1 16
./graphGenerator D6_part1.txt 20000000 400000000 1 17
./graphGenerator D6_part2.txt 20000000 400000000 1 18
./graphGenerator D6_part3.txt 20000000 400000000 1 19
./graphGenerator D7_part1.txt 20000000 500000000 1 20
./graphGenerator D7_part2.txt 20000000 600000000 1 21
./graphGenerator D8_part1.txt 30000000 600000000 1 22
./graphGenerator D8_part2.txt 20000000 500000000 1 23
./graphGenerator D9_part1.txt 10000000 400000000 1 24
./graphGenerator D9_part2.txt 10000000 500000000 1 25
./graphGenerator D10_part1.txt 10000000 300000000 1 26
./graphGenerator D10_part2.txt 10000000 300000000 1 27

./mergegraph D1_part1.txt D1.txt 1 0 1
./mergegraph D1_part2.txt D1.txt 0 100000000 1
./mergegraph D1_part3.txt D1.txt 0 200000000 1
./mergegraph D1_part4.txt D1.txt 0 300000000 1

./mergegraph D1.txt D2.txt 1 0 0
./mergegraph D2_part1.txt D2.txt 0 460000000 1
./mergegraph D2_part2.txt D2.txt 0 500000000 1
./mergegraph D2_part3.txt D2.txt 0 540000000

./mergegraph D2.txt D3.txt 1 0 0
./mergegraph D3_part1.txt D3.txt 0 630000000 1
./mergegraph D3_part2.txt D3.txt 0 650000000 1
./mergegraph D3_part3.txt D3.txt 0 670000000 1

./mergegraph D3.txt D4.txt 1 0 0
./mergegraph D4_part1.txt D4.txt 0 700000000 1
./mergegraph D4_part2.txt D4.txt 0 720000000 1
./mergegraph D4_part3.txt D4.txt 0 740000000 1

./mergegraph D4.txt D5.txt 1 0 0
./mergegraph D5_part1.txt D5.txt 0 770000000 1
./mergegraph D5_part2.txt D5.txt 0 790000000 1
./mergegraph D5_part3.txt D5.txt 0 810000000 1

./mergegraph D5.txt D6.txt 1 0 0
./mergegraph D6_part1.txt D6.txt 0 820000000 1
./mergegraph D6_part2.txt D6.txt 0 840000000 1
./mergegraph D6_part3.txt D6.txt 0 860000000 1

./mergegraph D6.txt D7.txt 1 0 0
./mergegraph D7_part1.txt D7.txt 0 880000000 1
./mergegraph D7_part2.txt D7.txt 0 900000000 1

./mergegraph D7.txt D8.txt 1 0 0
./mergegraph D8_part1.txt D8.txt 0 920000000 1
./mergegraph D8_part2.txt D8.txt 0 950000000 1

./mergegraph D8.txt D9.txt 1 0 0
./mergegraph D9_part1.txt D9.txt 0 970000000 1
./mergegraph D9_part2.txt D9.txt 0 980000000 1

./mergegraph D9.txt D10.txt 1 0 0
./mergegraph D10_part1.txt D10.txt 0 990000000 1
./mergegraph D10_part2.txt D10.txt 0 1000000000 1
