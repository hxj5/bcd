# subclone_inference.py - subclonal structure inference.

# Inputs
# * CalicoST - spot-wise clone label.
# * CopyKAT - hclustering results;
# * InferCNV - cell x gene CNA expression matrix;
# * Numbat - cell-wise clone label;
# * XClone - cell-wise clone label;


import os
from bcd.subclonal_structure import (
    subclonal_structure_main,
    CalicoST, CopyKAT, InferCNV, Numbat, XClone
)


sample_id = 'HCC3N_simu'
root_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n2t_aleloh/base'
in_dir = os.path.join(root_dir, 'cna_calling')
out_dir = os.path.join(root_dir, 'subclone_inference')

tool_list = []

clone_label_fn = os.path.join(in_dir, 'calicost/clone3_rectangle0_w1.0/clone_labels.tsv')
if os.path.exists(clone_label_fn):
    calicost = CalicoST(
        clone_label_fn = clone_label_fn
    )
    tool_list.append(calicost)

hclust_fn = os.path.join(in_dir, 'copykat/%s_copykat_clustering_results.rds' % sample_id)
if os.path.exists(hclust_fn):
    copykat = CopyKAT(
        hclust_fn = hclust_fn
    )
    tool_list.append(copykat)

obj_fn = os.path.join(in_dir, 'infercnv/run.final.infercnv_obj')
if os.path.exists(obj_fn):
    infercnv = InferCNV(
        obj_fn = obj_fn
    )
    tool_list.append(infercnv)

clone_post_fn = os.path.join(in_dir, 'numbat/cna/clone_post_2.tsv')
if os.path.exists(clone_post_fn):
    numbat = Numbat(
        clone_post_fn = clone_post_fn
    )
    tool_list.append(numbat)

clone_post_fn = os.path.join(in_dir, 'xclone/xclone_cna/data/HCC3N_simu_cell_clone_posteriors_modified-from-vb.tsv')
if os.path.exists(clone_post_fn):
    xclone = XClone(
        clone_post_fn = clone_post_fn,
        run_id = "4states-cluster-joint"
    )
    tool_list.append(xclone)

clone_post_fn = os.path.join(in_dir, 'xclone_aleloh/xclone_cna/data/HCC3N_simu_cell_clone_posteriors_modified-from-vb.tsv')
if os.path.exists(clone_post_fn):
    xclone = XClone(
        clone_post_fn = clone_post_fn,
        run_id = "5states-cluster-joint"
    )
    tool_list.append(xclone)

clone_post_fn = os.path.join(in_dir, 'xclone_aleloh_cluster-baf/xclone_cna/data/HCC3N_simu_cell_clone_posteriors_modified-from-vb.tsv')
if os.path.exists(clone_post_fn):
    xclone = XClone(
        clone_post_fn = clone_post_fn,
        run_id = "5states-cluster-baf"
    )
    tool_list.append(xclone)
    
subclonal_structure_main(
    sid = sample_id,
    tool_list = tool_list,
    out_dir = out_dir,
    truth_fn = os.path.join(root_dir, 'gen_data/simu/4_rs/rs.cell_anno.tsv'),
    ref_labels = 'normal',
    n_cluster = 2,
    merge_clusters = False,
    verbose = True
)
