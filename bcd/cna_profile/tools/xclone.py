# xclone.py


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



class XClone(Tool):
    def __init__(self, combine_fn):
        super().__init__(
            tid = "XClone",
            has_gain = True,
            has_loss = True,
            has_loh = True
        )
        self.combine_fn = combine_fn
        
        
    def extract(
        self,
        out_fn_list, 
        cna_type_list,
        verbose = False
    ):
        """Extract XClone prob matrix and convert it to python object.

        Parameters
        ----------
        out_fn_list : list of str
            Output ".h5ad" files storing the cell x gene matrix, each per cna type.
        cna_type_list : list of str
            A list of CNA types, each in {"gain", "loss", "loh"}.
        verbose : bool, default False
            Whether to show detailed logging information.

        Returns
        -------
        Void.
        """
        return extract_cna_prob(
            combine_fn = self.combine_fn,
            out_fn_list = out_fn_list, 
            cna_type_list = cna_type_list,
            verbose = verbose
        )



def extract_cna_prob(
    combine_fn,
    out_fn_list, 
    cna_type_list,
    verbose = False
):
    if verbose:
        info("check args ...")

    assert_e(combine_fn)
    assert len(out_fn_list) > 0
    assert len(cna_type_list) == len(out_fn_list)
    
    for cna_type in cna_type_list:
        assert cna_type in ("gain", "loss", "loh")
    
    if verbose:
        info("load XClone object ...")

    xclone_adata = ad.read_h5ad(combine_fn)
    
    if verbose:
        info("XClone adata shape = %s." % str(xclone_adata.shape))

    prob = xclone_adata.layers['prob1_merge']
    # prob_merge = np.stack([copy_loss, loh, copy_neutral, copy_gain], axis = -1)

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
            mtx = prob[:, :, 0]  # copy_loss
        elif cna_type == "loh":
            mtx = prob[:, :, 1]  # loh
        elif cna_type == "gain":
            mtx = prob[:, :, 3]  # copy_gain
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
        
    return(out_fn_list)
