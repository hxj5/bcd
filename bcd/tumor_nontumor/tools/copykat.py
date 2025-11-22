# copykat.py


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
from logging import info, error
from logging import warning as warn
from .base import Tool
from ..utils.base import assert_e
from ..utils.io import save_h5ad



class CopyKAT(Tool):
    def __init__(self, ploidy_pred_fn):
        """CopyKAT object.
        
        Parameters
        ----------
        ploidy_pred_fn : str
            Path to TSV file with columns 'cell.names' and 'copykat.pred'.
        """
        super().__init__(tid = "CopyKAT")
        self.ploidy_pred_fn = ploidy_pred_fn

        
    def predict(self, out_fn):
        """Process CopyKAT predictions of tumor vs. non-tumor.
        
        Read CopyKAT TSV file, rename columns to 'barcode' and 'prediction', 
        convert 'diploid' to 'normal' and 'aneuploid' to 'tumor', 
        and save to a new TSV file.

        Saves a TSV file with columns: 
        - ``barcode``, ``prediction`` ('normal' or 'tumor').
        """
        return extract_tumor_prediction(
            ploidy_pred_fn = self.ploidy_pred_fn,
            out_fn = out_fn,
            delimiter = '\t'
        )
        
        
        
def extract_tumor_prediction(
    ploidy_pred_fn,
    out_fn,
    delimiter = '\t'
):
    # Check args and load data.
    if not os.path.exists(ploidy_pred_fn):
        raise ValueError(f"TSV file not found at {ploidy_pred_fn}")
        
    df = pd.read_csv(ploidy_pred_fn, delimiter = delimiter)

    required_cols = ['cell.names', 'copykat.pred']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")

    # sometimes value 'not.defined' exists in column `copykat.pred`.
    valid_preds = {'diploid', 'aneuploid'}
    if not set(df['copykat.pred']).issubset(valid_preds):
        invalid_preds = set(df['copykat.pred']) - valid_preds
        warn(f"Invalid prediction values found: {invalid_preds}."  \
             f"Expected: {valid_preds}")
        df = df.loc[df['copykat.pred'].isin(valid_preds), :].copy()
        warn("%d cells left after removing invalid predictions!" % \
             df.shape[0])


    # Create output DataFrame
    result_df = pd.DataFrame({
        'barcode': df['cell.names'],
        'prediction': df['copykat.pred'].replace(
            {'diploid': 'normal', 'aneuploid': 'tumor'}),
        'ploidy_pred': df['copykat.pred']
    })

    
    # Save to TSV
    result_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Processed predictions saved to '{out_fn}'.")


    # Print summary
    n_cells = len(df)
    n_tumor = sum(result_df['prediction'] == 'tumor')
    info(f"Processed {n_cells} cells.")
    info(f"Number of tumor cells: {n_tumor}")
    info(f"Number of normal cells: {n_cells - n_tumor}")
    
    return(out_fn)
