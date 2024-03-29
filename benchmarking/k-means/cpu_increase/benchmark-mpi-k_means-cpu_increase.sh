#!/bin/bash
set -e

array=(1 2 4 8 12 16 20 24 28 30)
for nproc in "${array[@]}"; do
    nnodes=1
    ncpus=$nproc # to see max resource utilization
    input_file="/home/eleonora.renz/hpc4ds-project/benchmarking/k-means/heavy/Credit_Data_Heavy.csv"
    filename=benchmark-mpi-k_means-cpu_increase-$nnodes-$ncpus-$nproc.sh
    echo "#!/bin/bash

#PBS -l select=$nnodes:ncpus=$ncpus:mem=2gb 

# set max execution time
#PBS -l walltime=0:10:00

#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n $nproc ./hpc4ds-project/mpi-k_means $nnodes $ncpus $input_file" >$filename

    chmod u+x $filename

    if (($nproc % 2 == 0)); then # --> queue limit is 30 per user
        echo "Sleeping to not exceed queue limit"
        sleep 30
    fi

    for runs in {1..5}; do # to average the values later on
        qsub $filename
    done
done
