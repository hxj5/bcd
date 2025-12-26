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
            'clone_opt'.
        """
        super().__init__(tid = "Numbat")
        self.clone_post_fn = clone_post_fn
        

    def predict(self, out_fn, verbose = False):
        """Predict subclonal structure from Numbat output.

        Saves a TSV file with columns: `barcode`, `prediction`
        to out_dir/numbat_predictions.tsv.
        """
        return extract_clone_labels(
            clone_post_fn = self.clone_post_fn,
            out_fn = out_fn,
            barcode_col = 'cell',
            clone_col = 'clone_opt',
            delimiter = '\t',
            verbose = verbose
        )



def extract_clone_labels(
    clone_post_fn,
    out_fn,
    barcode_col = 'cell',
    clone_col = 'clone_opt',
    delimiter = '\t',
    verbose = False
):
    # Check args.
    assert_e(clone_post_fn)

    df = pd.read_csv(clone_post_fn, delimiter = delimiter)
    required_cols = [barcode_col, clone_col]
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")


    # Filter out rows with empty clone labels.
    initial_n_cells = len(df)
    df = df.dropna(subset = [clone_col])
    n_cells = len(df)
    n_removed = initial_n_cells - n_cells
    if n_removed > 0:
        warn("Removed %d cells with empty '%s' values." %  \
            (n_removed, clone_col))
    if n_cells == 0:
        error("No cells remain after removing empty '%s' values." % clone_col)
        raise ValueError


    # Format labels.
    labels_old = np.unique(df[clone_col])
    label_new = np.arange(len(labels_old))
    label_map = {o:n for o, n in zip(labels_old, label_new)}


    # Save to TSV
    result_df = pd.DataFrame({
        'barcode': df[barcode_col].to_numpy(),
        'prediction': df[clone_col].map(label_map),
        'clone_label': df[clone_col]
    })
    result_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Predictions saved to {out_fn}.")


    # Print summary
    df = pd.read_csv(out_fn, sep = '\t')
    labels = np.unique(df['prediction'])
    n_cells = len(df)
    info(f"Processed {n_cells} cells after filtering.")
    info("In total %d clones." % len(labels))

    return(out_fn)
