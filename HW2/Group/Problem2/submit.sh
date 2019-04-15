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
{ time ./parallel 8; } >sample.out.txt 2>&1
