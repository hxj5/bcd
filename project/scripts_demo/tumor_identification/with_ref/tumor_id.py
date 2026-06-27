# tumor_id.py - tumor identification.

# Inputs
# * CalicoST - spot-wise tumor prop.
# * CopyKAT - label of cell ploidy ('diploid' vs. 'aneuploid');
# * InferCNV - cell x gene CNA expression matrix;
# * Numbat - cell-wise 'tumor' vs. 'normal' assignment;
# * XClone - cell-wise 'tumor' vs. 'normal' assignment;


import os
from bcd.tumor_nontumor import (
    tumor_nontumor_main,
    CalicoST, CopyKAT, InferCNV, Numbat, XClone
)


sample_id = 'HCC3N_simu'
root_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1r1n1t/vary_tumor_prop/700n300t'
in_dir = os.path.join(root_dir, 'cna_calling')
out_dir = os.path.join(root_dir, 'tumor_identification/with_ref')

tool_list = []

tumor_prop_fn = os.path.join(in_dir, 'calicost/clone3_rectangle0_w1.0/clone_labels.tsv')
if os.path.exists(tumor_prop_fn):
    calicost = CalicoST(
        tumor_prop_fn = tumor_prop_fn
    )
    tool_list.append(calicost)

ploidy_pred_fn = os.path.join(in_dir, 'copykat/%s_copykat_prediction.txt' % sample_id)
if os.path.exists(ploidy_pred_fn):
    copykat = CopyKAT(
        ploidy_pred_fn = ploidy_pred_fn
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

xclone_tumor_pred_fn = os.path.join(in_dir, 'xclone/xclone_cna/data/xclone_tumor_predictions.tsv')
if os.path.exists(xclone_tumor_pred_fn):
    xclone = XClone(
        xclone_tumor_pred_fn = xclone_tumor_pred_fn
    )
    tool_list.append(xclone)
    
tumor_nontumor_main(
    sid = sample_id,
    tool_list = tool_list,
    out_dir = out_dir,
    truth_fn = os.path.join(root_dir, 'gen_data/matrix/spot_anno.tsv'),
    tumor_labels = 'tumor',
    ref_labels = 'ref',
    verbose = True
)
