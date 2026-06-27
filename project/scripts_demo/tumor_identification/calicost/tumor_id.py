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
out_dir = os.path.join(root_dir, 'tumor_identification/calicost')

tool_list = []

tumor_prop_fn = os.path.join(in_dir, 'calicost/clone3_rectangle0_w1.0/clone_labels.tsv')
if os.path.exists(tumor_prop_fn):
    calicost_ref_purity = CalicoST(
        tumor_prop_fn = tumor_prop_fn,
        run_id = 'ref_purity'
    )
    tool_list.append(calicost_ref_purity)
    
tumor_prop_fn = os.path.join(in_dir, 'calicost_noref/clone3_rectangle0_w1.0/clone_labels.tsv')
if os.path.exists(tumor_prop_fn):
    calicost_noref_purity = CalicoST(
        tumor_prop_fn = tumor_prop_fn,
        run_id = 'noref_purity'
    )
    tool_list.append(calicost_noref_purity)
    
tumor_prop_fn = os.path.join(in_dir, 'calicost_skippurity/clone3_rectangle0_w1.0/clone_labels.tsv')
if os.path.exists(tumor_prop_fn):
    calicost_ref_nopurity = CalicoST(
        tumor_prop_fn = tumor_prop_fn,
        run_id = 'ref_nopurity'
    )
    tool_list.append(calicost_ref_nopurity)
    
tumor_prop_fn = os.path.join(in_dir, 'calicost_noref_skippurity/clone3_rectangle0_w1.0/clone_labels.tsv')
if os.path.exists(tumor_prop_fn):
    calicost_noref_nopurity = CalicoST(
        tumor_prop_fn = tumor_prop_fn,
        run_id = 'noref_nopurity'
    )
    tool_list.append(calicost_noref_nopurity)

    
tumor_nontumor_main(
    sid = sample_id,
    tool_list = tool_list,
    out_dir = out_dir,
    truth_fn = os.path.join(root_dir, 'gen_data/matrix/spot_anno.tsv'),
    tumor_labels = 'tumor',
    ref_labels = 'ref',
    verbose = True
)
