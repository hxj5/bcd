# xclone_rdr.py


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
from logging import info, error
from logging import warning as warn
from sklearn.cluster import KMeans
from .base import Tool
from ..utils.base import assert_e
from ..utils.io import save_h5ad



class XCloneRDR(Tool):
    def __init__(self, adata_fn):
        """XClone object.
        
        Parameters
        ----------
        adata_fn : str
            Path to XClone adata file.
        """
        super().__init__(tid = "XClone_RDR")
        self.adata_fn = adata_fn
        

    def predict(self, out_fn, random_state = 123, verbose = False):
        """Predict tumor cells from XClone output.
        
        UPDATE ...

        Read XClone TSV file, apply K-means clustering on `p_cnv` column 
        to classify cells as 'normal' or 'tumor', and save predictions to 
        a TSV file in out_dir.

        Parameters
        ----------
        random_state : int
            Random seed for K-means (default: 123).

        Saves a TSV file with columns: `barcode`, `prediction`
        ('normal' or 'tumor') to out_dir/numbat_predictions.tsv.
        """
        return predict_tumor_from_cnv_prob(
            cnv_prob_fn = self.cnv_prob_fn,
            out_fn = out_fn,
            barcode_col = 'cell',
            p_cnv_col = 'p_cnv',
            delimiter = '\t',
            random_state = random_state,
            verbose = verbose
        )



def predict_tumor_from_expression(
    adata_fn,
    out_fn,
    layer = None,
    ref_expr = 0,
    linkage_method = 'ward',
    linkage_metric = 'euclidean',
    fcluster_criterion = 'maxclust',
    cna_score_how = 'mad',
    verbose = False
):
    """Predict tumor cells from XClone RDR output.
    
    Read AnnData, perform tumor classification using hierarchical 
    clustering on adata RDR, and save into TSV file: 
    - columns `barcode`, `prediction`.
    
    Saves TSV files: xclone_rdr_predictions.tsv.
    """
    # Check args.
    assert_e(adata_fn)
    
    
    # Perform tumor classification
    adata = ad.read_h5ad(adata_fn)
    X = None
    if layer is None:
        X = adata.X.copy()
    else:
        X = adata.layers[layer].copy()
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



def __predict_tumor_from_expression(
    X, 
    ref_expr = 0,
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
