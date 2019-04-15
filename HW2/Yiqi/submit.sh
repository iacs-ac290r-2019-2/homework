#!/bin/bash

#SBATCH --job-name=test
##SBATCH --account=ac290r
##SBATCH --reservation=ac290r_gpu
#SBATCH --partition=fas_gpu
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-00:05
##SBATCH -o job_%j.out
##SBATCH -e job_%j.err

hostname
echo $CUDA_VISIBLE_DEVICES
{ time ./parallel 4; } >trial_4.out 2>&1
{ time ./parallel 8; } >trial_8.out 2>&1
{ time ./parallel 16; } >trial_16.out 2>&1
{ time ./parallel 32; } >trial_32.out 2>&1
{ time ./parallel 64; } >trial_64.out 2>&1
{ time ./parallel 128; } >trial_128.out 2>&1
