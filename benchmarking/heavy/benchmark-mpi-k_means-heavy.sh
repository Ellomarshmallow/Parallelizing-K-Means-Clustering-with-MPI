#!/bin/bash
set -e

array=(1 2 4 8 16 32 64 128)
for nproc in "${array[@]}"; do
    nnodes=1
    ncpus=2
    input_file="/home/eleonora.renz/hpc4ds-project/benchmarking/heavy/Credit_Data_Heavy.csv"
    filename=benchmark-mpi-k_means-heavy-$nnodes-$ncpus-$nproc.sh
    echo "#!/bin/bash

#PBS -l select=$nnodes:ncpus=$ncpus:mem=2gb 

# set max execution time
#PBS -l walltime=0:10:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n $nproc ./hpc4ds-project/mpi-k_means-elli $nnodes $ncpus $input_file" >$filename

    chmod u+x $filename

    for runs in {1..3}; do # to average the values later on --> queue limit is 30 per user
        qsub $filename
    done
done
