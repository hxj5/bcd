#!/bin/bash
#PBS -N cellsnp_mode2a
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=80g,walltime=30:00:00
#PBS -o cellsnp_mode2a.out
#PBS -e cellsnp_mode2a.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi


# PUT YOUR CODE HERE
/usr/bin/time -v  cellsnp-lite  \
    -s  /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/sccnasim-cs_screadsim/align/cellranger/count/outs/possorted_genome_bam.bam   \
    -b  /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/sccnasim-cs_screadsim/simu/final/barcodes.tsv        \
    -O  /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/sccnasim-cs_screadsim/cellsnp_mode2a  \
    -p  20  \
    -f  /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/refseq/genome.fa  \
    --minCOUNT  11  \
    --minMAF  0.1    \
    --gzip

set +ux
conda deactivate
echo All Done!

