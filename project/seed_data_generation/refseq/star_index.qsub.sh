#!/bin/bash
#PBS -N star_index
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=100g,walltime=100:00:00
#PBS -o star_index.out
#PBS -e star_index.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
  work_dir=$PBS_O_WORKDIR
fi
work_dir=`cd $work_dir && cd .. && pwd`


# PUT YOUR CODE HERE
STAR  \
    --runThreadN  10  \
    --runMode  genomeGenerate  \
    --genomeDir  /groups/cgsd/xianjie/result/sccnasim/HCC3N_4289spot/data/star  \
    --genomeFastaFiles  /groups/cgsd/xianjie/data/refseq/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa  \
    --sjdbGTFfile  /groups/cgsd/xianjie/data/refseq/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf  \
    --sjdbOverhang  97


set +ux
conda deactivate
echo All Done!

