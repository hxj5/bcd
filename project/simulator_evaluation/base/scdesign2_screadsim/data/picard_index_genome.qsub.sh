#!/bin/bash
#PBS -N picard_index_genome
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=120g,walltime=20:00:00
#PBS -o picard_index_genome.out
#PBS -e picard_index_genome.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -eux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi


# PUT YOUR CODE HERE
cd /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/scdesign2_screadsim/data/refgenome

picard CreateSequenceDictionary \
    -R  genome.fa \
    -O  genome.fa.dict

set +ux
conda deactivate
echo All Done!

