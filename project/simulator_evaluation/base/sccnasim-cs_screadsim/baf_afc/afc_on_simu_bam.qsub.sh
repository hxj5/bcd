#!/bin/bash
#PBS -N afc_on_simu_bam
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=80g,walltime=15:00:00
#PBS -o afc_on_simu_bam.out
#PBS -e afc_on_simu_bam.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi


# PUT YOUR CODE HERE
/usr/bin/time -v  python  $work_dir/afc_on_simu_bam.py


set +ux
conda deactivate
echo All Done!
