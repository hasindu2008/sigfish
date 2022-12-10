#!/bin/bash

g++ -O2 -o ./UCR_DTW ./UCR_DTW.cpp
./UCR_DTW ref.txt query.txt 250 0.1

gcc -O2 -o ./cdtw ./cdtw.c ./dtw_main.c -lm
./cdtw ref.txt query.txt 29898 250