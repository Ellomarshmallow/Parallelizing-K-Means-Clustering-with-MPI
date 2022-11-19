#!/bin/bash

#PBS -l select=2:ncpus=2:mem=2gb

# set max execution time
#PBS -l walltime=1:00:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 4 ./hpc4ds-project/mpi-k_means-elli
