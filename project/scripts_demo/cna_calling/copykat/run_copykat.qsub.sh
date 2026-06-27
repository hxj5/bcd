#!/bin/bash
#PBS -N run_copykat
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=120g,walltime=20:00:00
#PBS -o run_copykat.out
#PBS -e run_copykat.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi


# PUT YOUR CODE HERE
/usr/bin/time -v  python $work_dir/run_copykat.py

set +ux
conda deactivate
echo All Done!

