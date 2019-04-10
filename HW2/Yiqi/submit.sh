#!/bin/bash

#SBATCH --job-name=test
##SBATCH --account=ac290r
##SBATCH --reservation=ac290r_gpu
#SBATCH --partition=fas_gpu
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-00:01
#SBATCH -o job_%j.out
#SBATCH -e job_%j.err

# Can run commands here, or use ‘echo’ to look at CUDA_VISIBLE_DEVICES
hostname
echo $CUDA_VISIBLE_DEVICES
time ./parallel 8
