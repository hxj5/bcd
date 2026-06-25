#!/bin/bash
#PBS -N run_cellranger
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=120g,walltime=30:00:00
#PBS -o run_cellranger.out
#PBS -e run_cellranger.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi


# PUT YOUR CODE HERE
/usr/bin/time -v  /home/xianjie/tools/cellranger-4.0.0/bin/cellranger count  \
    --id=HCC3N_simu        \
    --transcriptome=/groups/cgsd/xianjie/data/refseq/refdata-cellranger-GRCh38-3.0.0   \
    --fastqs=/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/sccnasim-cs_screadsim/align/cellranger/data \
    --sample=HCC3Nsimu     \
    --localcores=5         \
    --localmem=100         \
    --force-cells=1200

set +ux
conda deactivate
echo All Done!

