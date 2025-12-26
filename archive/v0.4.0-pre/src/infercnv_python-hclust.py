# infercnv.py


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
from logging import info, error
from logging import warning as warn
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.cluster import KMeans
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
        ref_expr = 1,
        linkage_method = 'ward',
        linkage_metric = 'euclidean',
        fcluster_criterion = 'maxclust',
        cna_score_how = 'mad',
        verbose = False
    ):
        os.makedirs(out_dir, exist_ok = True)
        out_fn = os.path.join(out_dir, "%s_predictions.tsv" % self.tid.lower())
        res = predict_tumor_from_expression(
            obj_fn = self.obj_fn,
            out_fn = out_fn,
            tmp_dir = os.path.join(out_dir, "r2py"),
            ref_expr = ref_expr,
            linkage_method = linkage_method,
            linkage_metric = linkage_metric,
            fcluster_criterion = fcluster_criterion,
            cna_score_how = cna_score_how,
            verbose = verbose
        )
        return out_fn



def predict_tumor_from_expression(
    obj_fn,
    out_fn,
    tmp_dir,
    ref_expr = 1,
    linkage_method = 'ward',
    linkage_metric = 'euclidean',
    fcluster_criterion = 'maxclust',
    cna_score_how = 'mad',
    verbose = False
):
    """Predict tumor cells from inferCNV output.
    
    Read AnnData, perform tumor classification using hierarchical 
    clustering on adata.X, and save into TSV file: 
    - columns `barcode`, `prediction`.
    
    Saves TSV files: infercnv_predictions.tsv.
    """
    # Check args.
    assert_e(obj_fn)
    os.makedirs(tmp_dir, exist_ok = True)
    
    
    # convert rds to adata
    adata_fn = os.path.join(tmp_dir, 'r2py.h5ad')
    extract_cna_expression(obj_fn, adata_fn, tmp_dir = tmp_dir)
    
    
    # Perform tumor classification
    adata = ad.read_h5ad(adata_fn)
    X = adata.X.copy()
    tumor_pred, scores = __predict_tumor_from_expression(
        X = X, 
        ref_expr = ref_expr,
        linkage_method = linkage_method,
        linkage_metric = linkage_metric,
        fcluster_criterion = fcluster_criterion,
        cna_score_how = cna_score_how
    )
    
    
    # Save predictions TSV
    df = pd.DataFrame({
        'barcode': adata.obs['cell'],
        'prediction': tumor_pred,
        'cna_score': scores
    })
    df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Predictions saved to {out_fn}.")
    
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
    
    
    
def __predict_tumor_from_expression(
    X, 
    ref_expr = 1,
    linkage_method = 'ward',
    linkage_metric = 'euclidean',
    fcluster_criterion = 'maxclust',
    cna_score_how = 'mad'
):
    """Cluster cells based on expression matrix using hierarchical clustering.
    
    Cluster cells based on expression matrix using hierarchical clustering 
    with (by default) Ward linkage and Euclidean distance, forcing 2 clusters.
    Label the cluster with the lowest aberration score
    (mean absolute deviation from median expression) as 'normal' and 
    the other as 'tumor'.
    
    Parameters
    ----------
    X : matrix-like
        The cell x gene expression matrix.
    ref_expr : float
        The reference expression value for calculating the deviation.
        Typically use the median expression value.
    
    Returns
    -------
    np.array
        The predicted labels: 'normal' or 'tumor'.
    """
    n_clusters = 2
    n_cells, n_genes = X.shape
    
    # Apply hierarchical clustering with Ward linkage and Euclidean distance
    Z = linkage(X, method = linkage_method, metric = linkage_metric)
    cluster_labels = fcluster(Z, t = n_clusters, criterion = fcluster_criterion)
    cluster_labels = cluster_labels - 1      # To make labels: 0 or 1
    
    
    # Assign cluster labels based on aberration score.
    # - Compute aberration score as mean absolute deviation from
    #   median expression per gene.
    # ref_expr = np.median(X, axis = 0)      # Median per gene
    scores = None
    if cna_score_how == 'mad':
        scores = np.mean(np.abs(X - ref_expr), axis = 1)   # Shape: (n_cells,)
    elif cna_score_how == 'md':
        scores = np.mean(X - ref_expr, axis = 1)   # Shape: (n_cells,)
    else:
        raise ValueError("unknown cna_score_how '%s'!" % cna_score_how)
    mean_scores = [np.mean(scores[cluster_labels == i]) for i in [0, 1]]
    normal_cluster = np.argmin(mean_scores)
    tumor_pred = np.where(cluster_labels == normal_cluster, 'normal', 'tumor')
    
    
    # Print summary
    info(f"Processed {n_cells} cells and {n_genes} genes.")
    info(f"Mean aberration scores per cluster: {mean_scores}")
    info(f"Cluster labels: {['normal', 'tumor']}")
    for label in np.unique(tumor_pred):
        count = np.sum(tumor_pred == label)
    info(f"Number of cells in {label}: {count}")
    
    return (tumor_pred, scores)
