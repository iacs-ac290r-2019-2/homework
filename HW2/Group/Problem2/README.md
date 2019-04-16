# Instruction

Upload this whole directory to your HPC. 

To compile, run
```
$ module load cuda
$ make
```

To submit, run
```
sbatch submit.sh
```

The default runs a 8x8 grid with `MAX_THREADS_DIM`=18. The numeric errors and the elapsed time will be reported in `sample.out.txt`. 
