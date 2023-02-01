#!/bin/bash
#SBATCH --account=rrg-yjliu
#SBATCH --ntasks=256                      # number of MPI processes
#SBATCH --mem-per-cpu=2G                  # memory; default unit is megabytes
#SBATCH --time=0-11:05                    # time (DD-HH:MM)
srun ../../bin/exe_solver                 # mpirun or mpiexec also work
