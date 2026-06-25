#!/bin/bash
#PBS -N cellsnp_mode2a
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=60g,walltime=20:00:00
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
    -s  /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/seed/xf_tag/bam/seed.xf_tag.possort.bam   \
    -b  /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/cell_anno/barcodes.tsv        \
    -O  /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/seed/xf_tag/cellsnp_mode2a  \
    -p  20  \
    -f  /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/refseq/genome.fa  \
    --minCOUNT  11  \
    --minMAF  0.1    \
    --gzip

set +ux
conda deactivate
echo All Done!

