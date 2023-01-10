#!/bin/bash
set -e

for nproc in {10..10..5}; do
    nnodes=1
    ncpus=2
    filename=benchmark-$nnodes-$ncpus-$nproc.sh
    echo "#!/bin/bash

#PBS -l select=$nnodes:ncpus=$ncpus:mem=2gb 

# set max execution time
#PBS -l walltime=1:00:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n $nproc ./hpc4ds-project/mpi-k_means-elli $nnodes $ncpus" >$filename

    chmod u+x $filename

    for runs in {1..10}; do # --> avg in python
        qsub $filename
    done
done
