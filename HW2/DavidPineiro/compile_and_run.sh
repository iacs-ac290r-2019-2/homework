#!/bin/bash

gcc -Wall -Wextra -O2 ./LBM_2D_channel.c -o LBM_2D_channel.x
time ./LBM_2D_channel.x > ux.dat

