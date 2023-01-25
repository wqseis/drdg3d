#!/bin/bash
#SBATCH --account=def-yjliu
#SBATCH --ntasks=64                     # number of MPI processes
#SBATCH --mem-per-cpu=256              # memory; default unit is megabytes
#SBATCH --time=0-01:05                 # time (DD-HH:MM)
srun ../../bin/exe_solver                 # mpirun or mpiexec also work
