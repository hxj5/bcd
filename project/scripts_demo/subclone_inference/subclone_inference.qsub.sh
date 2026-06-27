#!/bin/bash
#PBS -N subclone_inference
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=120g,walltime=10:00:00
#PBS -o subclone_inference.out
#PBS -e subclone_inference.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi


# PUT YOUR CODE HERE
/usr/bin/time -v  python $work_dir/subclone_inference.py

set +ux
conda deactivate
echo All Done!

