#!/bin/bash
# -*- coding: utf-8 -*-

./graphGenerator D1_part1.bin 100000000 400000000 0 1
./graphGenerator D1_part2.bin 100000000 400000000 0 2
./graphGenerator D1_part3.bin 100000000 400000000 0 3
./graphGenerator D1_part4.bin 160000000 500000000 0 4
./graphGenerator D2_part1.bin 50000000 400000000 0 5
./graphGenerator D2_part2.bin 50000000 400000000 0 6
./graphGenerator D2_part3.bin 70000000 500000000 0 7
./graphGenerator D3_part1.bin 20000000 400000000 0 8
./graphGenerator D3_part2.bin 20000000 400000000 0 9
./graphGenerator D3_part3.bin 30000000 500000000 0 10
./graphGenerator D4_part1.bin 20000000 400000000 0 11
./graphGenerator D4_part2.bin 20000000 400000000 0 12
./graphGenerator D4_part3.bin 30000000 400000000 0 13
./graphGenerator D5_part1.bin 20000000 400000000 0 14
./graphGenerator D5_part2.bin 20000000 400000000 0 15
./graphGenerator D5_part3.bin 10000000 400000000 0 16
./graphGenerator D6_part1.bin 20000000 400000000 0 17
./graphGenerator D6_part2.bin 20000000 400000000 0 18
./graphGenerator D6_part3.bin 20000000 400000000 0 19
./graphGenerator D7_part1.bin 20000000 500000000 0 20
./graphGenerator D7_part2.bin 20000000 600000000 0 21
./graphGenerator D8_part1.bin 30000000 600000000 0 22
./graphGenerator D8_part2.bin 20000000 500000000 0 23
./graphGenerator D9_part1.bin 10000000 400000000 0 24
./graphGenerator D9_part2.bin 10000000 500000000 0 25
./graphGenerator D10_part1.bin 10000000 300000000 0 26
./graphGenerator D10_part2.bin 10000000 300000000 0 27

./normalize ./D1_part1.bin
./normalize ./D1_part2.bin
./normalize ./D1_part3.bin
./normalize ./D1_part4.bin
./normalize ./D2_part1.bin
./normalize ./D2_part2.bin
./normalize ./D2_part3.bin
./normalize ./D3_part1.bin
./normalize ./D3_part2.bin
./normalize ./D3_part3.bin
./normalize ./D4_part1.bin
./normalize ./D4_part2.bin
./normalize ./D4_part3.bin
./normalize ./D5_part1.bin
./normalize ./D5_part2.bin
./normalize ./D5_part3.bin
./normalize ./D6_part1.bin
./normalize ./D6_part2.bin
./normalize ./D6_part3.bin
./normalize ./D7_part1.bin
./normalize ./D7_part2.bin
./normalize ./D8_part1.bin
./normalize ./D8_part2.bin
./normalize ./D9_part1.bin
./normalize ./D9_part2.bin
./normalize ./D10_part1.bin
./normalize ./D10_part2.bin
