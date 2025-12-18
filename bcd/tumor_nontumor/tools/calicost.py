# calicost.py


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



class CalicoST(Tool):
    def __init__(self, tumor_prop_fn):
        """CalicoST Object.
        
        Parameters
        ----------
        tumor_prop_fn : str
            Path to CalicoST TSV file containing columns: 
            ``BARCODES``, ``clone_label``, and ``tumor_proportion``.
        """
        super().__init__(tid = "CalicoST")
        self.tumor_prop_fn = tumor_prop_fn


    def predict(self, out_fn, verbose = False):
        """Predict tumor cells from CalicoST output.
        
        Predict tumor cells from CalicoST output using K-means on 
        `tumor_proportion`, excluding cells with empty tumor_proportion,
        and save to a TSV file with columns 'barcode' and 'prediction' 
        ('tumor' or 'normal').

        Saves a TSV file with columns: `barcode`, `prediction` ('tumor' or 
        'normal'), and `tumor_proportion` for cells with valid 
        tumor_proportion.
        """
        return predict_tumor_from_prop(
            tumor_prop_fn = self.tumor_prop_fn,
            out_fn = out_fn,
            prop_col = 'tumor_proportion',
            delimiter = '\t',
            random_state = 123,
            verbose = verbose
        )



def predict_tumor_from_prop(
    tumor_prop_fn,
    out_fn,
    prop_col,
    delimiter = '\t',
    random_state = 123,
    verbose = False
):
    """
    Parameters
    ----------
    random_state : int
        Random seed for K-means (default: 123).
    """
    n_clusters = 2       # tumor and normal (non-tumor).

    # Check args and load data.
    df = pd.read_csv(tumor_prop_fn, delimiter = delimiter)

    required_cols = ['BARCODES', 'clone_label', prop_col]
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")


    # Filter out rows with empty tumor_proportion
    initial_n_cells = len(df)
    df = df.dropna(subset = [prop_col])
    n_cells = len(df)
    n_removed = initial_n_cells - n_cells
    if n_removed > 0:
        warn("Removed %d cells with empty '%s' values." %  \
            (n_removed, prop_col))
    if n_cells == 0:
        error("No cells remain after removing empty '%s' values." % prop_col)
        raise ValueError

        
    # Apply K-means clustering
    tumor_proportions = df[prop_col].values.reshape(-1, 1)
    kmeans = KMeans(n_clusters = n_clusters, random_state = random_state)
    cluster_labels = kmeans.fit_predict(tumor_proportions)

    
    # Identify tumor cluster (higher mean tumor_proportion)
    mean_scores = np.zeros(n_clusters)
    for cluster in range(n_clusters):
        cluster_cells = cluster_labels == cluster
        mean_scores[cluster] = np.mean(tumor_proportions[cluster_cells])
    tumor_cluster = np.argmax(mean_scores)
    predictions = np.where(
        cluster_labels == tumor_cluster, 'tumor', 'normal')


    # Save to TSV
    result_df = pd.DataFrame({
        'barcode': df['BARCODES'],
        'prediction': predictions,
        'tumor_proportion': tumor_proportions
    })
    result_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Predictions saved to '{out_fn}'.")

    
    # Print summary
    
    # below codes only work for 2 clusters
    # Compute threshold (midpoint between cluster centers)
    cluster_centers = kmeans.cluster_centers_.flatten()
    low_center, high_center = sorted(cluster_centers)
    threshold = (low_center + high_center) / 2

    n_tumor = np.sum(predictions == 'tumor')
    info(f"Processed {n_cells} cells after filtering.")
    info(f"%s tumor_proportion cluster centers: %s." % \
         (self.tid, cluster_centers))
    info(f"Selected threshold: {threshold:.4f} " \
         "(cells > threshold classified as tumor)")
    info(f"Number of tumor cells: {n_tumor}.")
    info(f"Number of non-tumor cells: {n_cells - n_tumor}.")
    
    return(out_fn)
