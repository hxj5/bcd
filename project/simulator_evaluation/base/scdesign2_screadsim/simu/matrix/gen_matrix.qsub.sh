#!/bin/bash
#PBS -N gen_matrix
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=60g,walltime=20:00:00
#PBS -o gen_matrix.out
#PBS -e gen_matrix.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -ux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi


# PUT YOUR CODE HERE
/usr/bin/time -v  python $work_dir/gen_matrix.py

set +ux
conda deactivate
echo All Done!

