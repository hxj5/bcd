# numbat.py


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



class Numbat(Tool):
    def __init__(self, cnv_prob_fn):
        """Numbat object.
        
        Parameters
        ----------
        cnv_prob_fn : str
            Path to Numbat TSV file with columns including 'cell' and 'p_cnv'.
        """
        super().__init__(tid = "Numbat")
        self.cnv_prob_fn = cnv_prob_fn
        

    def predict(self, out_fn, random_state = 123):
        """Predict tumor cells from Numbat output.

        Read Numbat TSV file, apply K-means clustering on `p_cnv` column 
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
            random_state = random_state
        )



def predict_tumor_from_cnv_prob(
    cnv_prob_fn,
    out_fn,
    barcode_col = 'cell',
    p_cnv_col = 'p_cnv',
    delimiter = '\t',
    random_state = 123
):
    n_clusters = 2

    # Check args.
    assert_e(cnv_prob_fn)

    df = pd.read_csv(cnv_prob_fn, delimiter = delimiter)
    required_cols = [barcode_col, p_cnv_col]
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")


    # Filter out rows with empty p_cnv
    initial_n_cells = len(df)
    df = df.dropna(subset = [p_cnv_col])
    n_cells = len(df)
    n_removed = initial_n_cells - n_cells
    if n_removed > 0:
        warn(f"Removed {n_removed} cells with empty {p_cnv_col} values.")
    if n_cells == 0:
        error(f"No cells remain after removing empty '%s' values." % \
              p_cnv_col)
        raise ValueError

    
    # Apply K-means clustering
    p_cnv_values = df[p_cnv_col].values.reshape(-1, 1)
    kmeans = KMeans(n_clusters = n_clusters, random_state = random_state)
    cluster_labels = kmeans.fit_predict(p_cnv_values)


    # Identify tumor cluster (higher mean p_cnv)
    mean_scores = np.zeros(n_clusters)
    for cluster in range(n_clusters):
        cluster_cells = cluster_labels == cluster
        mean_scores[cluster] = np.mean(p_cnv_values[cluster_cells])
    tumor_cluster = np.argmax(mean_scores)
    predictions = np.where(
        cluster_labels == tumor_cluster, 'tumor', 'normal')


    # Save to TSV
    result_df = pd.DataFrame({
        'barcode': df[barcode_col],
        'prediction': predictions,
        'p_cnv': p_cnv_values
    })
    result_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Predictions saved to {out_fn}.")


    # Print summary

    # CHECK ME: below codes only work for 2 clusters.
    # Compute threshold (midpoint between cluster centers)
    cluster_centers = kmeans.cluster_centers_.flatten()
    low_center, high_center = sorted(cluster_centers)
    threshold = (low_center + high_center) / 2

    n_tumor = np.sum(predictions == 'tumor')
    info(f"Processed {n_cells} cells after filtering.")
    info(f"Numbat p_cnv cluster centers: {cluster_centers}")
    info(f"Selected threshold: {threshold:.4f} " \
         "(cells > threshold classified as tumor)")
    info(f"Number of tumor cells: {n_tumor}")
    info(f"Number of normal cells: {n_cells - n_tumor}")

    return(out_fn)
