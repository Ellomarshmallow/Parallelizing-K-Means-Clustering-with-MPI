#!/bin/bash

#PBS -l select=2:ncpus=2:mem=2gb

# set max execution time
#PBS -l walltime=0:10:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 2 ./hpc4ds-project/mpi-k_means 2 2