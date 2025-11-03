# xclone_rdr.py


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



class XCloneRDR(Tool):
    def __init__(self, obj_path, rdr_fn = None):
        """
        rdr_fn : str
            File storing the read depth ratio (RDR) matrix.
            If None, do not use RDR.
        """
        super().__init__(
            tid = "XCloneRDR",
            obj_path = obj_path,
            has_gain = True,
            has_loss = True,
            has_loh = False            
        )
        self.rdr_fn = rdr_fn

        
    def extract(
        self,
        out_fn_list, 
        cna_type_list,
        tmp_dir,
        verbose = False
    ):
        """Extract XClone RDR matrix and convert it to python object.

        Parameters
        ----------
        out_fn_list : list of str
            Output ".h5ad" files storing the cell x gene matrix, each per cna type.
        cna_type_list : list of str
            A list of CNA types, each in {"gain", "loss", "loh"}.
        tmp_dir : str
            The folder to store temporary data.
        verbose : bool, default False
            Whether to show detailed logging information.

        Returns
        -------
        Void.
        """
        obj_fn, rdr_fn = self.obj_path, self.rdr_fn
        
        if verbose:
            info("check args ...")

        assert_e(obj_fn)
        assert_e(rdr_fn)

        assert len(out_fn_list) > 0
        assert len(cna_type_list) == len(out_fn_list)

        for cna_type in cna_type_list:
            assert cna_type in ("gain", "loss")

        os.makedirs(tmp_dir, exist_ok = True)

        if verbose:
            info("load XClone object ...")

        xclone_adata = ad.read_h5ad(obj_fn)

        if verbose:
            info("XClone adata shape = %s." % str(xclone_adata.shape))

        rdr_adata = ad.read_h5ad(rdr_fn)

        if verbose:
            info("RDR adata shape = %s." % str(rdr_adata.shape))

        # Prepare obs and var DataFrames
        cells = xclone_adata.obs.index.tolist()
        genes = xclone_adata.var['GeneName'].tolist()

        obs_df = pd.DataFrame(data = dict(cell = cells))
        var_df = pd.DataFrame(data = dict(gene = genes))

        # iterate over cna types.
        for cna_type, out_fn in zip(cna_type_list, out_fn_list):
            if verbose:
                info("process cna_type '%s' ..." % cna_type)

            if cna_type == "loss":
                mtx = xclone_adata.layers['posterior_mtx'][:, :, 0]  # copy_loss
            elif cna_type == "gain":
                mtx = xclone_adata.layers['posterior_mtx'][:, :, 1]  # copy_gain
            else:
                raise ValueError(f"Error: unknown cnv type '{cna_type}'.")

            adata = ad.AnnData(
                X = mtx,
                obs = obs_df,
                var = var_df
            )
            save_h5ad(adata, out_fn)
            
            if verbose:
                info("saved adata shape = %s." % str(adata.shape))

            del adata
            gc.collect()
