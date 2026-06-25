# afc_on_simu_bam.py - run `afc` on the simulated BAM file.


import os
import sys
from sccnasim.afc.main import afc_wrapper
from sccnasim.main import add_cell_anno


if len(sys.argv) < 3:
    print("Usage: %s <run_dir> <out_dir>" % sys.argv[0])
    sys.exit(1)
    
run_dir = sys.argv[1]
out_dir = sys.argv[2]

cell_anno_fn = os.path.join(run_dir, "simu/4_rs/rs.cell_anno.tsv")
assert os.path.exists(cell_anno_fn)


afc_wrapper(
    sam_fn = os.path.join(run_dir, "simu/4_rs/bam/rs.possorted.bam"), 
    barcode_fn = os.path.join(run_dir, "simu/4_rs/rs.samples.tsv"), 
    feature_fn = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/gene_anno/gene_anno.tsv", 
    phased_snp_fn = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/snp/HCC3N.phased.vcf.gz", 
    out_dir = out_dir,
    ncores = 10,
    min_count = 1, 
    min_maf = 0,
    strandness = "forward"
)

count_fn = os.path.join(out_dir, "afc.counts.h5ad")

count_fn_new = count_fn.replace(".h5ad", ".cell_anno.h5ad")
add_cell_anno(
    count_fn = count_fn,
    cell_anno_fn = cell_anno_fn,
    out_count_fn = count_fn_new
)
