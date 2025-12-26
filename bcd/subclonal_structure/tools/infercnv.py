# infercnv.py


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
from logging import info, error
from logging import warning as warn
from .base import Tool
from ..utils.base import assert_e, exe_cmdline
from ..utils.io import save_h5ad



class InferCNV(Tool):
    def __init__(self, obj_fn):
        """InferCNV object.
        
        Parameters
        ----------
        obj_fn : str
            File storing the inferCNV object. 
            Typically using the "MCMC_inferCNV_obj.rds".
        """
        super().__init__(tid = "inferCNV")
        self.obj_fn = obj_fn
        
        
    def predict(
        self,
        out_dir,
        k,
        dist = 'euclidean',
        hclust = 'ward.D2',
        verbose = False
    ):
        os.makedirs(out_dir, exist_ok = True)
        out_fn = os.path.join(out_dir, "%s_predictions.tsv" % self.tid.lower())
        res = predict_subclones_from_expression(
            obj_fn = self.obj_fn,
            out_fn = out_fn,
            k = k,
            tmp_dir = os.path.join(out_dir, "r2py"),
            dist = dist,
            hclust = hclust,
            verbose = verbose
        )
        return out_fn



def predict_subclones_from_expression(
    obj_fn,
    out_fn,
    k,
    tmp_dir,
    dist = 'euclidean',
    hclust = 'ward.D2',
    verbose = False
):
    # Check args.
    assert_e(obj_fn)
    os.makedirs(tmp_dir, exist_ok = True)
    
    
    # convert rds to adata
    adata_fn = os.path.join(tmp_dir, 'r2py.h5ad')
    extract_cna_expression(obj_fn, adata_fn, tmp_dir = tmp_dir)
    
    
    # Perform subclone identification.
    s  = ""
    s += '''# Identify subclones from expression.\n'''
    s += '''\n'''
    s += '''obj <- readRDS("%s")\n''' % obj_fn
    s += '''mtx <- obj@expr.data\n'''
    s += '''mtx <- t(mtx)         # cell x gene matrix\n'''
    s += '''hc <- hclust(dist(mtx, method = "%s"), method = "%s")\n''' % (dist, hclust)  
    s += '''label <- cutree(tree = hc, k = %d)\n''' % k
    s += '''df <- data.frame(\n'''
    s += '''    barcode = gsub(".", "-", names(label), fixed = TRUE),\n'''
    s += '''    prediction = label - 1\n'''
    s += ''')\n'''
    s += '''write.table(\n'''
    s += '''    df,\n'''
    s += '''    file = "%s",\n''' % out_fn
    s += '''    sep = "\\t",\n'''
    s += '''    row.names = FALSE,\n'''
    s += '''    col.names = TRUE\n'''
    s += ''')\n'''
    s += '''\n'''
    
    script_fn = os.path.join(tmp_dir, "predict_subclones_from_expression.R")
    with open(script_fn, "w") as fp:
        fp.write(s)


    # run the R script.
    if verbose:
        info("run the R script to predict subclones from expression ...")
    exe_cmdline("Rscript %s" % script_fn)
    
    
    # Save to TSV
    info(f"Processed predictions saved to '{out_fn}'.")


    # Print summary
    df = pd.read_csv(out_fn, sep = '\t')
    labels = np.unique(df['prediction'])
    n_cells = len(df)
    info(f"Processed {n_cells} cells after filtering.")
    info("In total %d clones." % len(labels))
    
    return(out_fn)



def extract_cna_expression(obj_fn, out_fn, tmp_dir, verbose = False):
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
    str
        Converted adata file.
    """
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
    
    s  = ""
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
        
    return(out_fn)
