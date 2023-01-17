#!/bin/bash

#PBS -l select=3:ncpus=3:mem=12gb
#PBS -l walltime=00:15:00
#PBS -q short_cpuQ
#PBS -M telucero@gmail.com



module load mpich-3.2

# Run the MPI program with 4 processes

for np in {1..5}
do
mpirun.actual -np 30 ./HPC4DS/go 3 3 /home/taylor.lucero/HPC4DS/Credit_Data_Heavy.csv 
done

