#!/bin/bash

gcc -Wall -Wextra -O2 ./LBM_2D_channel_v1.c -o LBM_2D_channel_v1.x
time ./LBM_2D_channel_v1.x > ux.dat

