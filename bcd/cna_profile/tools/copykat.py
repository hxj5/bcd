# copykat.py


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
import scipy as sp
from logging import info
from .base import Tool
from ..utils.base import assert_e, exe_cmdline
from ..utils.gscale import reg2gene
from ..utils.io import load_gene_anno, save_h5ad



class CopyKAT(Tool):
    def __init__(self, expr_mtx_fn):
        super().__init__(
            tid = "CopyKat",
            has_gain = True,
            has_loss = True,
            has_loh = False
        )
        self.expr_mtx_fn = expr_mtx_fn


    def extract(
        self,
        out_fn,
        verbose = False
    ):
        """
        Extract CNV data from CopyKAT output and save as AnnData object.

        Parameters
        ----------
        out_fn : str
            Output .h5ad file to save the matrix. 
        verbose : bool, default False
            Whether to print verbose output.

        Returns
        -------
        dict
            {"mtx": cell x gene matrix (numpy array), "overlap": None}
        """
        return extract_cna_expression(
            expr_mtx_fn = self.expr_mtx_fn,
            out_fn = out_fn,
            verbose = verbose
        )



def extract_cna_expression(
    expr_mtx_fn,
    out_fn,
    verbose = False
):
    if verbose:
        info(f"Reading CopyKAT file: {expr_mtx_fn}")
    assert_e(expr_mtx_fn)

    mtx = pd.read_csv(expr_mtx_fn, sep = "\t", header = 0, dtype = str)
    mtx.index = mtx["hgnc_symbol"]
    mtx = mtx.iloc[:, 7:]              # Remove first 7 columns.
    mtx = mtx.T                        # cell x gene

    if verbose:
        info(f"CopyKAT matrix shape: {mtx.shape}")

    adata = ad.AnnData(
        X = mtx.values.astype(float),
        obs = pd.DataFrame(data = dict(cell = mtx.index)),
        var = pd.DataFrame(data = dict(gene = mtx.columns))
    )
    save_h5ad(adata, out_fn)
    
    if verbose:
        info(f"Saved AnnData to {out_fn}")
        
    return(out_fn)
