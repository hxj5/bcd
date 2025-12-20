# numbat.py


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
from logging import info, error
from logging import warning as warn
from .base import Tool
from ..utils.base import assert_e



class Numbat(Tool):
    def __init__(self, clone_post_fn):
        """Numbat object.
        
        Parameters
        ----------
        clone_post_fn : str
            Path to Numbat TSV file with columns including 'cell' and 
            'compartment_opt'.
        """
        super().__init__(tid = "Numbat")
        self.clone_post_fn = clone_post_fn
        

    def predict(self, out_fn, verbose = False):
        """Predict tumor cells from Numbat output.

        Saves a TSV file with columns: `barcode`, `prediction`
        ('normal' or 'tumor') to out_dir/numbat_predictions.tsv.
        """
        return extract_tumor(
            clone_post_fn = self.clone_post_fn,
            out_fn = out_fn,
            barcode_col = 'cell',
            label_col = 'compartment_opt',
            delimiter = '\t',
            verbose = verbose
        )



def extract_tumor(
    clone_post_fn,
    out_fn,
    barcode_col = 'cell',
    label_col = 'compartment_opt',
    p_cnv_col = 'p_cnv',
    delimiter = '\t',
    verbose = False
):
    # Check args.
    assert_e(clone_post_fn)

    df = pd.read_csv(clone_post_fn, delimiter = delimiter)
    required_cols = [barcode_col, label_col, p_cnv_col]
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")


    # Filter out rows with empty or invalid labels.
    n_cells_init = len(df)
    df = df.dropna(subset = [label_col])
    df = df.loc[df[label_col].isin(['normal', 'tumor'])].copy()
    n_cells = len(df)
    n_removed = n_cells_init - n_cells
    if n_removed > 0:
        warn(f"Removed %d cells with empty/invalid '%s' values." % \
            (n_removed, label_col))
    if n_cells == 0:
        error(f"No cells remain after removing empty/invalid '%s' values." % \
              label_col)
        raise ValueError


    # Save to TSV
    result_df = pd.DataFrame({
        'barcode': df[barcode_col].to_numpy(),
        'prediction': df[label_col].to_numpy(),
        'p_cnv': df[p_cnv_col].to_numpy()
    })
    result_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Predictions saved to {out_fn}.")
    
    
    # Print Summary.
    df = pd.read_csv(out_fn, sep = '\t')
    n_cells = df.shape[0]
    n_tumor = np.sum(df['prediction'] == 'tumor')
    n_normal = np.sum(df['prediction'] == 'normal')
    info("Number of all_cells=%d; tumor_cells=%d; normal_cells=%d." % \
        (n_cells, n_tumor, n_normal))

    return(out_fn)
