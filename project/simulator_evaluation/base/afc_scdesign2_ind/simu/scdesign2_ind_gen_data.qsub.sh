#!/bin/bash
#PBS -N scdesign2_ind_gen_data
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=120g,walltime=20:00:00
#PBS -o scdesign2_ind_gen_data.out
#PBS -e scdesign2_ind_gen_data.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi


# PUT YOUR CODE HERE
/usr/bin/time -v  Rscript $work_dir/scdesign2_ind_gen_data.R

set +ux
conda deactivate
echo All Done!

