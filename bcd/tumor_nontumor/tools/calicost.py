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
from ..utils.io import save_h5ad



class CalicoST(Tool):
    def __init__(
        self,
        obj_path, 
        out_dir,
        proportion_col = 'tumor_proportion', 
        delimiter = '\t'
    ):
        """Initialize CalicoST tool with directory containing input files.
        
        Parameters
        ----------
        obj_path : str
            Path to CalicoST TSV file containing columns: 
            ``BARCODES``, ``clone_label``, and ``tumor_proportion``.
        proportion_col : str, default 'tumor_proportion'
            Column name in `obj_path` for tumor proportions.
        delimiter : str, default: '\t'
            Delimiter for TSV file.
        """
        super().__init__(
            tid = "CalicoST",
            obj_path = obj_path,
            out_dir = out_dir
        )
        self.proportion_col = proportion_col
        self.delimiter = delimiter

        
    def predict(
        self,
        n_clusters = 2,
        random_state = 42,
        verbose = False
    ):
        """Predict tumor cells from CalicoST output.
        
        Predict tumor cells from CalicoST output using K-means on 
        `tumor_proportion`, excluding cells with empty tumor_proportion,
        and save to a TSV file with columns 'barcode' and 'prediction' 
        ('tumor' or 'normal').

        Parameters:
        -----------
        n_clusters : int
            Number of clusters for K-means (default: 2 for tumor/non-tumor).
        random_state : int
            Random seed for K-means (default: 42).

        Returns:
        --------
        None
            Saves a TSV file with columns: barcode, prediction ('tumor' or 
            'normal') for cells with valid tumor_proportion.
        """
        tsv_path = self.obj_path
        out_dir = self.out_dir
        proportion_col = self.proportion_col
        delimiter = self.delimiter


        # Check args.
        if not out_dir:
            raise ValueError("out_dir must be provided to save predictions.")

        df = pd.read_csv(tsv_path, delimiter = delimiter)

        required_cols = ['BARCODES', 'clone_label', proportion_col]
        if not all(col in df.columns for col in required_cols):
            raise ValueError(f"TSV must contain columns: {required_cols}")


        # Filter out rows with empty tumor_proportion
        initial_n_cells = len(df)
        df = df.dropna(subset=[proportion_col])
        n_cells = len(df)
        n_removed = initial_n_cells - n_cells
        if n_removed > 0:
            warn("Removed %d cells with empty '%s' values." %  \
                (n_removed, proportion_col))
        if n_cells == 0:
            error("No cells remain after removing empty '%s' values." % \
                  proportion_col)
            raise ValueError

            
        # Apply K-means clustering
        tumor_proportions = df[proportion_col].values.reshape(-1, 1)
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

        
        # Compute threshold (midpoint between cluster centers)
        cluster_centers = kmeans.cluster_centers_.flatten()
        low_center, high_center = sorted(cluster_centers)
        threshold = (low_center + high_center) / 2

        result_df = pd.DataFrame({
            'barcode': df['BARCODES'],
            'prediction': predictions
        })

        
        # Save to TSV
        fn = os.path.join(
            out_dir, '%s_predictions.tsv' % self.tid.lower())
        result_df.to_csv(predictions_path, sep = '\t', index = False)
        info(f"Predictions saved to '{fn}'.")

        
        # Print summary
        n_tumor = np.sum(predictions == 'tumor')
        info(f"Processed {n_cells} cells after filtering.")
        info(f"CalicoST tumor_proportion cluster centers: {cluster_centers}")
        info(f"Selected threshold: {threshold:.4f} "   \
             "(cells > threshold classified as tumor)")
        info(f"Number of tumor cells: {n_tumor}")
        info(f"Number of non-tumor cells: {n_cells - n_tumor}")
