#!/bin/bash

make
gcc -Wall -O2 -g -I include -I slow5lib/include examples/ex.c -o ex lib/libsigfish.a slow5lib/lib/libslow5.a -lz -lm -lpthread
