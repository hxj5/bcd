import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
from logging import info, warning
from .base import Tool
from ..utils.base import assert_e
from ..utils.io import save_h5ad

class CalicoST(Tool):
    def __init__(self, obj_dir):
        """Initialize CalicoST tool with directory containing input files.
        
        Parameters
        ----------
        obj_dir : str
            Directory containing 'cnv_genelevel.tsv' and 'clone_labels.tsv' files.
        """
        super().__init__(
            tid="CalicoST",
            obj_path=obj_dir,
            has_gain=True,
            has_loss=True,
            has_loh=True
        )
        self.cnv_file = os.path.join(obj_dir, "cnv_genelevel.tsv")
        self.clone_file = os.path.join(obj_dir, "clone_labels.tsv")

    def extract(
        self,
        out_fn_list,
        cna_type_list,
        tmp_dir,
        verbose=False
    ):
        """Extract CalicoST data and convert it to cell x gene probability matrices.
        Probabilities are scaled by tumor proportion for each cell.
        
        Parameters
        ----------
        out_fn_list : list of str
            Output ".h5ad" files storing the cell x gene matrix, each per CNA type.
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
        if verbose:
            info("Checking arguments...")
        assert_e(self.cnv_file)
        assert_e(self.clone_file)
        assert len(out_fn_list) > 0
        assert len(cna_type_list) == len(out_fn_list)
        for cna_type in cna_type_list:
            assert cna_type in ("gain", "loss", "loh")

        os.makedirs(tmp_dir, exist_ok=True)

        if verbose:
            info("Loading CalicoST data...")
        # Load input files
        cnv_df = pd.read_csv(self.cnv_file, sep='\t')
        clone_df = pd.read_csv(self.clone_file, sep='\t')

        # Validate input files
        assert 'gene' in cnv_df.columns, "cnv_genelevel.tsv must contain 'gene' column"
        assert all(col in clone_df.columns for col in ['BARCODES', 'clone_label', 'tumor_proportion']), \
            "clone_labels.tsv must contain 'BARCODES', 'clone_label', and 'tumor_proportion' columns"

        # Handle duplicate barcodes and genes
        if clone_df['BARCODES'].duplicated().any():
            warning(f"Found {clone_df['BARCODES'].duplicated().sum()} duplicate barcodes in clone_labels.tsv. Keeping first occurrence.")
            clone_df = clone_df.drop_duplicates(subset='BARCODES', keep='first')
        if cnv_df['gene'].duplicated().any():
            warning(f"Found {cnv_df['gene'].duplicated().sum()} duplicate genes in cnv_genelevel.tsv. Keeping first occurrence.")
            cnv_df = cnv_df.drop_duplicates(subset='gene', keep='first')

        # Extract genes, cells, and tumor proportions
        genes = cnv_df['gene'].tolist()
        cells = clone_df['BARCODES'].tolist()
        clone_labels = clone_df['clone_label'].values
        tumor_proportions = clone_df['tumor_proportion'].values
        unique_clones = np.unique(clone_labels)

        # Initialize obs and var DataFrames
        obs_df = pd.DataFrame(data=dict(cell=cells))
        var_df = pd.DataFrame(data=dict(gene=genes))

        # Ensure unique indices
        # obs_df = obs_df.set_index('cell', verify_integrity=True)
        # var_df = var_df.set_index('gene', verify_integrity=True)

        # Process each CNA type
        for cna_type, out_fn in zip(cna_type_list, out_fn_list):
            if verbose:
                info(f"Processing CNA type '{cna_type}'...")

            # Initialize cell x gene matrix
            mtx = np.zeros((len(cells), len(genes)), dtype=np.float32)

            # Process each clone
            for clone in unique_clones:
                # Get cells belonging to this clone
                cell_mask = clone_df['clone_label'] == clone
                cell_indices = np.where(cell_mask)[0]

                # Get copy number data for this clone
                col_a = f'clone{clone} A'
                col_b = f'clone{clone} B'
                if col_a not in cnv_df.columns or col_b not in cnv_df.columns:
                    if verbose:
                        info(f"Skipping clone {clone}: missing A or B copy number columns")
                    continue

                # Compute total copy number (A + B)
                total_cn = cnv_df[col_a] + cnv_df[col_b]

                # Assign base probabilities based on CNA type
                if cna_type == "loss":
                    prob = (total_cn < 2).astype(np.float32)  # Loss: A+B < 2
                elif cna_type == "loh":
                    prob = ((total_cn == 2) & ((cnv_df[col_a] == 0) | (cnv_df[col_b] == 0))).astype(np.float32)  # LOH: A+B = 2 and (A=0 or B=0)
                elif cna_type == "gain":
                    prob = (total_cn > 2).astype(np.float32)  # Gain: A+B > 2
                else:
                    raise ValueError(f"Error: unknown CNA type '{cna_type}'.")

                # Scale probabilities by tumor proportion for each cell in the clone
                for cell_idx in cell_indices:
                    mtx[cell_idx, :] = prob * tumor_proportions[cell_idx]

            # Create AnnData object
            adata = ad.AnnData(
                X=mtx,
                obs=obs_df,
                var=var_df
            )

            # Save to h5ad file
            save_h5ad(adata, out_fn)
            if verbose:
                info(f"Saved adata shape = {adata.shape} for CNA type '{cna_type}'.")

            # Clean up
            del adata
            gc.collect()