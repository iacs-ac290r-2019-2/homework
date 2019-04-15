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

The numeric errors are listed in the stdout, and the timing report is in the stderr. For reference, you can find my results in `sample.out.txt` and `sample.err.txt`. 
