#!/bin/bash
#PBS -N subset_bam
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=100g,walltime=100:00:00
#PBS -o subset_bam.out
#PBS -e subset_bam.err

source ~/.bashrc
conda activate F

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
  work_dir=$PBS_O_WORKDIR
fi


# PUT YOUR CODE HERE
echo "Subset BAM ..."

sample=HCC3N_600spot
samtools view -h -b \
  -o  $work_dir/${sample}.possort.bam  \
  -D  CB:/groups/cgsd/xianjie/projects/cna-benchmark/dev/HCC3N_600spot/data/cell_anno/barcodes.tsv  \
  -@  5   \
  /groups/cgsd/xianjie/projects/cna-benchmark/dev/HCC3N_4289spot/data/bam/HCC3N.bam


echo "Index BAM ..."
samtools index \
  -@  5   \
  $work_dir/${sample}.possort.bam


set +ux
conda deactivate
echo All Done!
