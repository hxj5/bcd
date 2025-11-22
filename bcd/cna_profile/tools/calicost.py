# calicost.py


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
from logging import info, warning
from .base import Tool
from ..utils.base import assert_e
from ..utils.io import save_h5ad



class CalicoST(Tool):
    def __init__(self, cnv_fn, clone_fn):
        """Initialize CalicoST."""
        super().__init__(
            tid = "CalicoST",
            has_gain = True,
            has_loss = True,
            has_loh = True
        )
        self.cnv_fn = cnv_fn
        self.clone_fn = clone_fn

        
    def extract(self, out_fn_list, cna_type_list, verbose = False):
        """Extract CalicoST data and convert it to cell x gene probability 
        matrices.
        
        Probabilities are scaled by tumor proportion for each cell.
        
        Parameters
        ----------
        out_fn_list : list of str
            Output ".h5ad" files storing the cell x gene matrix, each per 
            CNA type.
        cna_type_list : list of str
            A list of CNA types, each in {"gain", "loss", "loh"}.
        verbose : bool, default False
            Whether to show detailed logging information.
        
        Returns
        -------
        Void.
        """
        return extract_cna_prob(
            cnv_fn = self.cnv_fn,
            clone_fn = self.clone_fn,
            out_fn_list = out_fn_list,
            cna_type_list = cna_type_list,
            verbose = verbose            
        )



def extract_cna_prob(
    cnv_fn,
    clone_fn,
    out_fn_list,
    cna_type_list,
    verbose = False
):
    if verbose:
        info("Checking arguments...")

    assert_e(cnv_fn)
    assert_e(clone_fn)
    assert len(out_fn_list) > 0
    assert len(cna_type_list) == len(out_fn_list)
    for cna_type in cna_type_list:
        assert cna_type in ("gain", "loss", "loh")

        
    # Load input data.
    if verbose:
        info("Loading CalicoST data ...")

    cnv_df = pd.read_csv(cnv_fn, sep = '\t')
    clone_df = pd.read_csv(clone_fn, sep = '\t')

    assert 'gene' in cnv_df.columns, \
        "cnv_genelevel.tsv must contain 'gene' column!"
    for col in ['BARCODES', 'clone_label', 'tumor_proportion']:
        assert col in clone_df.columns, \
            "clone_labels.tsv must contain '%s' column!" % col

        
    # Handle duplicate barcodes and genes
    if clone_df['BARCODES'].duplicated().any():
        n_dup = clone_df['BARCODES'].duplicated().sum()
        warning(f"Found {n_dup} duplicate barcodes in clone_labels.tsv. " \
                "Keeping first occurrence.")
        clone_df = clone_df.drop_duplicates(
            subset = 'BARCODES', keep = 'first')

    if cnv_df['gene'].duplicated().any():
        n_dup = cnv_df['gene'].duplicated().sum()
        warning(f"Found {n_dup} duplicate genes in cnv_genelevel.tsv. "  \
                "Keeping first occurrence.")
        cnv_df = cnv_df.drop_duplicates(
            subset = 'gene', keep = 'first')


    # Init output data.
    genes = cnv_df['gene'].tolist()
    cells = clone_df['BARCODES'].tolist()
    clone_labels = clone_df['clone_label'].values
    tumor_proportions = clone_df['tumor_proportion'].values
    unique_clones = np.unique(clone_labels)

    
    obs_df = pd.DataFrame(data = dict(cell = cells))
    var_df = pd.DataFrame(data = dict(gene = genes))

    # Ensure unique indices
    # obs_df = obs_df.set_index('cell', verify_integrity = True)
    # var_df = var_df.set_index('gene', verify_integrity = True)

    
    # Process each CNA type
    for cna_type, out_fn in zip(cna_type_list, out_fn_list):
        if verbose:
            info(f"Processing CNA type '{cna_type}'...")

        mtx = np.zeros((len(cells), len(genes)), dtype = np.float32)
        for clone in unique_clones:
            cell_mask = clone_df['clone_label'] == clone
            cell_indices = np.where(cell_mask)[0]

            col_a = f'clone{clone} A'
            col_b = f'clone{clone} B'
            assert col_a in cnv_df.columns, \
                'clone{clone}: missing A copy number in cnv_df!'
            assert col_b in cnv_df.columns, \
                'clone{clone}: missing B copy number in cnv_df!'

            total_cn = cnv_df[col_a] + cnv_df[col_b]

            # Assign base probabilities based on CNA type
            if cna_type == "loss":         # Loss: A+B < 2
                prob = (total_cn < 2)
            elif cna_type == "loh":        # LOH: A+B = 2 and (A=0 or B=0)
                prob = ((total_cn == 2) & \
                        ((cnv_df[col_a] == 0) | (cnv_df[col_b] == 0)))
            elif cna_type == "gain":       # Gain: A+B > 2
                prob = (total_cn > 2)  
            else:
                raise ValueError(f"Error: unknown CNA type '{cna_type}'.")
            prob = prob.astype(np.float32)

            # Scale probabilities by tumor proportion for each cell 
            # in the clone
            # CHECK ME: P(spot belongs to specific clone) == tumor_prop?
            for cell_idx in cell_indices:
                mtx[cell_idx, :] = prob * tumor_proportions[cell_idx]

        adata = ad.AnnData(
            X = mtx,
            obs = obs_df,
            var = var_df
        )

        save_h5ad(adata, out_fn)
        if verbose:
            info(f"Saved adata shape = '%s' for CNA type '%s'." %  \
                 (str(adata.shape), cna_type))

        del adata
        gc.collect()

    return(out_fn_list)
