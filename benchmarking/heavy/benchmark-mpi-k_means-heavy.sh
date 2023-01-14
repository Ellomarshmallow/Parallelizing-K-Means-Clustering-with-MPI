#!/bin/bash
set -e

array=(1 2 4 8 16 32 64 128)
for nproc in "${array[@]}"; do
    nnodes=1
    ncpus=2
    #dataset_size=10000
    filename=benchmark-heavy-$nnodes-$ncpus-$nproc.sh
    echo "#!/bin/bash

#PBS -l select=$nnodes:ncpus=$ncpus:mem=2gb 

# set max execution time
#PBS -l walltime=1:00:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n $nproc ./hpc4ds-project/mpi-k_means-elli $nnodes $ncpus" >$filename #missing $dataset_size because it didn't work yet

    chmod u+x $filename

    for runs in {1..3}; do # to average the values later on --> queue limit is 30 per user
        qsub $filename
    done
done
