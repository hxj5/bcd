#!/bin/bash
#PBS -N subset_bam
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=120g,walltime=20:00:00
#PBS -o subset_bam.out
#PBS -e subset_bam.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi


# PUT YOUR CODE HERE
samtools view -h -b -s 0.5 -@ 10 \
    -o  /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/downsample_coverage/ds_50perc/gen_data/subset_bam/possorted.bam  \
    /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/base/data/possorted.bam

samtools index -b \
    /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/downsample_coverage/ds_50perc/gen_data/subset_bam/possorted.bam


set +ux
conda deactivate
echo All Done!

