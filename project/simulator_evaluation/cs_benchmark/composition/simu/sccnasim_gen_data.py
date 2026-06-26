# sccnasim_gen_data.py


import os
from sccnasim import main_wrapper


dataset_data_dir = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data"
study_data_dir = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/cs_benchmark/composition/mixed/data"
out_dir = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/cs_benchmark/composition/mixed/simu"


main_wrapper(
    sam_fn = os.path.join(dataset_data_dir, "bam/HCC3N_600spot.possort.bam"),
    cell_anno_fn = os.path.join(dataset_data_dir, "cell_anno/spot_anno.tsv"),
    feature_fn = os.path.join(dataset_data_dir, "gene_anno/gene_anno.tsv"),
    phased_snp_fn = os.path.join(dataset_data_dir, "snp/HCC3N.phased.vcf.gz"),
    clone_anno_fn = os.path.join(study_data_dir, "clone_anno.tsv"), 
    cna_profile_fn = os.path.join(study_data_dir, "cna_profile.tsv"), 
    refseq_fn = os.path.join(dataset_data_dir, "refseq/genome.fa"), 
    out_dir = out_dir,
    libsize_ratio = 1.0,
    umi_len = 12,
    barcode_whitelist_fn = os.path.join(dataset_data_dir, "refapp/cellranger/737K-august-2016.txt"),
    ncores = 10,
    strandness = "forward"
)
