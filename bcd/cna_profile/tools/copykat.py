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
    def __init__(self, obj_path):
        super().__init__(
            tid = "CopyKat",
            obj_path = obj_path,
            has_gain = True,
            has_loss = True,
            has_loh = False
        )


    def extract(
        self,
        out_fn,
        tmp_dir,
        verbose = False
    ):
        """
        Extract CNV data from CopyKAT output and save as AnnData object.

        Parameters
        ----------
        out_fn : str
            Output .h5ad file to save the matrix.
        tmp_dir : str
            The folder to store temporary data. 
        verbose : bool, default False
            Whether to print verbose output.

        Returns
        -------
        dict
            {"mtx": cell x gene matrix (numpy array), "overlap": None}
        """
        obj_fn = self.obj_path
        
        if verbose:
            info(f"Reading CopyKAT file: {obj_fn}")
        assert_e(obj_fn)

        os.makedirs(tmp_dir, exist_ok = True)

        mtx = pd.read_csv(obj_fn, sep = "\t", header = 0, dtype = str)
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
