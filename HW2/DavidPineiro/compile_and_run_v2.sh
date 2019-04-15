#!/bin/bash

gcc -Wall -Wextra -O2 ./LBM_2D_channel_v2.c -o LBM_2D_channel_v2.x
time ./LBM_2D_channel_v2.x > ux.dat

