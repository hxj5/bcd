# truth.py - format input ground truth.


import anndata as ad
import functools
import gc
import numpy as np
import os
import pandas as pd
from logging import info, error
from ..utils.base import assert_e



def run_truth(
    truth_fn,
    out_fn,
    tumor_labels,
    verbose = True
):
    """Format ground truth.
    
    Read a TSV file containing ground truth annotations, extract barcode 
    and true label columns, and save to a TSV file.
    
    Saves a TSV file with columns: `barcode`, `annotation`.

    Parameters:
    -----------
    truth_fn : str
        The 2-column header-free TSV file with ground truth annotations:
        `barcode` and `annotation`.
    out_fn : str
        File to save the formatted ground truth.
    tumor_labels : str or list of str
        The cell type labels for tumor cells.
    """
    barcode_col = 'barcode'
    anno_col = 'annotation'

    # Check args.
    if not os.path.exists(truth_fn):
        raise ValueError(f"TSV file not found at {truth_fn}!")
    
    if isinstance(tumor_labels, str):
        tumor_labels = [tumor_labels]

    df = pd.read_csv(truth_fn, delimiter = '\t', header = None)
    df.columns = [barcode_col, anno_col]


    # Create output DataFrame
    out_df = pd.DataFrame({
        barcode_col: df[barcode_col]
    })
    out_df[anno_col] = 'normal'
    out_df.loc[df[anno_col].isin(tumor_labels), anno_col] = 'tumor'


    # Save to TSV
    out_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Ground truth saved to {out_fn}.")


    # Print summary
    n_cells = len(out_df)
    n_unique_labels = len(out_df[anno_col].unique())
    info(f"Extracted {n_cells} cells with ground truth annotations.")
    info(f"Unique labels in {anno_col}: {out_df[anno_col].unique()}.")
        
        
    res = dict(
        # out_fn : str
        #   Tumor label ground truth, TSV file.
        out_fn = out_fn
    )
    return(res)
