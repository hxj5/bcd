#!/bin/bash
#PBS -N subset_xf_tag
#PBS -q cgsd
#PBS -l nodes=1:ppn=3,mem=10g,walltime=10:00:00
#PBS -o subset_xf_tag.out
#PBS -e subset_xf_tag.err

source ~/.bashrc
conda activate SCSC

# run `set` after `source` & `conda activate` as the source file has an unbound variable
set -ux  

work_dir=`cd $(dirname $0) && pwd`
if [ -n "$PBS_O_WORKDIR" ]; then
    work_dir=$PBS_O_WORKDIR
fi


# PUT YOUR CODE HERE
out_bam=/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/sccnasim-cs_screadsim/align/cellranger/xf_tag/bam/sccnasim-cs_screadsim.xf_tag.possort.bam

/usr/bin/time -v  python $work_dir/subset_xf_tag.py  \
    /groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/sccnasim-cs_screadsim/align/cellranger/count/outs/possorted_genome_bam.bam   \
    $out_bam

samtools index  -@ 5  $out_bam


set +ux
conda deactivate
echo All Done!
