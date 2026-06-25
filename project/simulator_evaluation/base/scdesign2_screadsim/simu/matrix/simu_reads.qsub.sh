#!/bin/bash
#PBS -N simu_reads
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=200g,walltime=30:00:00
#PBS -o simu_reads.out
#PBS -e simu_reads.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -ux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi


# PUT YOUR CODE HERE
/usr/bin/time -v  python $work_dir/simu_reads.py  &>$work_dir/simu_reads.log

set +ux
conda deactivate
echo All Done!

