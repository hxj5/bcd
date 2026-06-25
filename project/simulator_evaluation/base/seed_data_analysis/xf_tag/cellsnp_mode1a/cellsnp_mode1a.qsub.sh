#!/bin/bash
#PBS -N cellsnp_mode1a
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=60g,walltime=20:00:00
#PBS -o cellsnp_mode1a.out
#PBS -e cellsnp_mode1a.err

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
    -R  /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/snp/HCC3N.phased.vcf.gz  \
    -O  /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/seed/xf_tag/cellsnp_mode1a  \
    -p  20  \
    --minCOUNT  0  \
    --minMAF  0    \
    --gzip

set +ux
conda deactivate
echo All Done!

