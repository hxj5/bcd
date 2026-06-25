#!/bin/bash
#PBS -N pp
#PBS -q cgsd
#PBS -l nodes=1:ppn=5,mem=60g,walltime=20:00:00
#PBS -o pp.out
#PBS -e pp.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -ux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi


# PUT YOUR CODE HERE
bam_fn=/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/sccnasim-cs_screadsim/data/HCC3N_600spot.possort.bam
out_dir=/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/sccnasim-cs_screadsim/simu/pp

mkdir $out_dir
samtools view $bam_fn -H > $out_dir/unprocess.header.sam

# create a bam file with the barcode embedded into the read name
time(cat <( cat $out_dir/unprocess.header.sam ) \
 <( samtools view $bam_fn | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
 | samtools view -bS - > $out_dir/processed.bam)

set +ux
conda deactivate
echo All Done!

