# cna_profile.py - CNA profile detection.

# Inputs
# * CalicoST - clonal gene copy number, spot-wise clone label and tumor prop.
# * CopyKAT - cell x gene CNA expression matrix;
# * InferCNV - cell x gene CNA expression matrix;
# * Numbat - cell x segment CNA probability;
# * XClone - cell x gene CNA probability;


import os
from bcd.cna_profile import (
    cna_profile_main,
    CalicoST, CopyKAT, InferCNV, Numbat, XClone
)


sample_id = 'HCC3N_simu'
root_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/downsample_coverage/ds_50perc'
in_dir = os.path.join(root_dir, 'cna_calling')
out_dir = os.path.join(root_dir, 'cna_prediction')

calicost = CalicoST(
    cnv_fn = os.path.join(in_dir, 'calicost/clone3_rectangle0_w1.0/cnv_genelevel.tsv'),
    clone_fn = os.path.join(in_dir, 'calicost/clone3_rectangle0_w1.0/clone_labels.tsv')
)
copykat = CopyKAT(
    expr_mtx_fn = os.path.join(in_dir, 'copykat/%s_copykat_CNA_raw_results_gene_by_cell.txt' % sample_id)
)
infercnv = InferCNV(
    obj_fn = os.path.join(in_dir, 'infercnv/run.final.infercnv_obj')
)
numbat = Numbat(
    joint_post_fn = os.path.join(in_dir, 'numbat/cna/joint_post_2.tsv')
)
xclone = XClone(
    combine_fn = os.path.join(in_dir, 'xclone/xclone_cna/data/combined_final.h5ad')
)

cna_profile_main(
    sid = sample_id,
    tool_list = [calicost, copykat, infercnv, numbat, xclone],
    out_dir = out_dir,
    truth_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/base/data/cna_truth.tsv',
    cell_anno_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/base/data/matrix/spot_anno.tsv',
    gene_anno_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/gene_anno/gene_anno.4column.tsv',
    ref_labels = ["ref"],
    cna_type_list = None,
    overlap_how = "isec-cells",
    max_n_cutoff = 1000,
    verbose = True
)
