# infercnv.py


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



        
class InferCNV(Tool):
    def __init__(self, obj_path):
        """
        obj_path : str
            File storing the inferCNV object. Typically using the
            "BayesNetOutput.HMMi6.hmm_mode-samples/MCMC_inferCNV_obj.rds".
        """
        super().__init__(
            tid = "inferCNV", 
            obj_path = obj_path,
            has_gain = True,
            has_loss = True,
            has_loh = False
        )


    def extract(self, out_fn, tmp_dir, verbose = False):
        """Extract inferCNV expression matrix and convert it to python object.

        Parameters
        ----------
        out_fn : str
            Output ".h5ad" file storing the cell x gene matrix.
        tmp_dir : str
            The folder to store temporary data.
        verbose : bool, default False
            Whether to show detailed logging information.

        Returns
        -------
        Void.
        """
        obj_fn = self.obj_path
        
        # check args.
        if verbose:
            info("check args ...")

        assert_e(obj_fn)
        os.makedirs(tmp_dir, exist_ok = True)


        # generate R scripts.
        if verbose:
            info("generate R scripts ...")

        cell_fn = os.path.join(tmp_dir, "barcodes.tsv")
        gene_fn = os.path.join(tmp_dir, "genes.tsv")
        mtx_fn = os.path.join(tmp_dir, "matrix.mtx")

        s = ""
        s += '''# extract cell x gene expression matrix from infercnv output.\n'''
        s += '''\n'''
        s += '''obj <- readRDS("%s")\n''' % obj_fn
        s += '''mtx <- obj@expr.data\n'''
        s += '''mtx <- t(mtx)         # cell x gene matrix\n'''
        s += '''\n'''  
        s += '''write(\n'''
        s += '''    rownames(mtx),\n'''
        s += '''    file = "%s"\n''' % cell_fn
        s += ''')\n'''
        s += '''\n'''
        s += '''write(\n'''
        s += '''    colnames(mtx),\n'''
        s += '''    file = "%s"\n''' % gene_fn
        s += ''')\n'''
        s += '''\n'''
        s += '''write.table(\n'''
        s += '''    mtx,\n'''
        s += '''    file = "%s",\n''' % mtx_fn
        s += '''    row.names = FALSE,\n'''
        s += '''    col.names = FALSE\n'''
        s += ''')\n'''
        s += '''\n'''

        script_fn = os.path.join(tmp_dir, "extract_infercnv.R")
        with open(script_fn, "w") as fp:
            fp.write(s)


        # run the R script to save expression matrix into file.
        if verbose:
            info("run the R script to save expression matrix into file ...")

        exe_cmdline("Rscript %s" % script_fn)


        # load matrix into anndata.
        if verbose:
            info("load matrix into anndata and save into .h5ad file ...")

        barcodes = pd.read_csv(cell_fn, header = None)
        barcodes.columns = ["cell"]

        genes = pd.read_csv(gene_fn, header = None)
        genes.columns = ["gene"]

        mtx = np.loadtxt(mtx_fn)

        adata = ad.AnnData(
            X = mtx,
            obs = barcodes,
            var = genes
        )
        save_h5ad(adata, out_fn)

        if verbose:
            info("saved adata shape = %s." % str(adata.shape))
