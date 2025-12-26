# calicost.py


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
from logging import info, error
from logging import warning as warn
from .base import Tool
from ..utils.base import assert_e



class CalicoST(Tool):
    def __init__(self, clone_label_fn):
        """CalicoST Object.
        
        Parameters
        ----------
        clone_label_fn : str
            Path to CalicoST TSV file containing columns: 
            ``BARCODES`` and ``clone_label``.
        """
        super().__init__(tid = "CalicoST")
        self.clone_label_fn = clone_label_fn


    def predict(self, out_fn, verbose = False):
        """Extract the subclonal labels from CalicoST output.

        Saves a TSV file with columns: `barcode`, `prediction`,
        and `clone_label`.
        """
        return extract_clonal_labels(
            clone_label_fn = self.clone_label_fn,
            out_fn = out_fn,
            delimiter = '\t',
            verbose = verbose
        )



def extract_clonal_labels(
    clone_label_fn,
    out_fn,
    delimiter = '\t',
    verbose = False
):
    # Check args and load data.
    df = pd.read_csv(clone_label_fn, delimiter = delimiter)

    required_cols = ['BARCODES', 'clone_label']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")


    # Filter out rows with empty clone labels.
    initial_n_cells = len(df)
    df = df.dropna(subset = ['clone_label'])
    n_cells = len(df)
    n_removed = initial_n_cells - n_cells
    if n_removed > 0:
        warn("Removed %d cells with empty 'clone_label' values." %  \
            (n_removed, ))
    if n_cells == 0:
        error("No cells remain after removing empty 'clone_label' values.")
        raise ValueError


    # Format labels.
    labels_old = np.unique(df['clone_label'])
    label_new = np.arange(len(labels_old))
    label_map = {o:n for o, n in zip(labels_old, label_new)}
    

    # Save to TSV
    result_df = pd.DataFrame({
        'barcode': df['BARCODES'],
        'prediction': df['clone_label'].map(label_map),
        'clone_label': df['clone_label']
    })
    result_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Predictions saved to '{out_fn}'.")

    
    # Print summary
    info(f"Processed {n_cells} cells after filtering.")
    info("In total %d clones." % len(labels_new))
    
    return(out_fn)
