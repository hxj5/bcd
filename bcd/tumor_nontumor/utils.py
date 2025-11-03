# Tumor non_tumor classification benchmark utils
# Author: Jiamu James Qiao
# Date: 2025-10

# tool specific processing:
# XClone: use the tumor_classify function, expression mtx -> prediction
# InferCNV: same function as xclone, need preprocessing
# Numbat: extract from R data, k means
# CopyKat: directly read output file
# CalicoST: k means on tumor proportion

# inputs:
# xclone: RDR / final adata
# infercnv: bcd output adata
# numbat


# xclone
import anndata as ad
import pandas as pd
import os
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

# truth
import pandas as pd
import os

def extract_ground_truth(
    tsv_path,
    out_dir='output',
    barcode_col='barcode',
    true_col='annotation',
    delimiter='\t'
):
    """
    Read a TSV file containing ground truth annotations, extract barcode and true label columns,
    and save to a TSV file in out_dir.

    Parameters:
    -----------
    tsv_path : str
        Path to TSV file with ground truth annotations.
    out_dir : str
        Directory to save the output TSV file (ground_truth.tsv).
    barcode_col : str
        Column name for cell barcodes (default: 'barcode').
    true_col : str
        Column name for true labels (default: 'annotation').
    delimiter : str
        Delimiter for the input TSV file (default: '\t').

    Returns:
    --------
    None
        Saves a TSV file with columns: barcode, annotation to out_dir/ground_truth.tsv.
    """
    # Validate input file
    if not os.path.exists(tsv_path):
        raise ValueError(f"TSV file not found at {tsv_path}")

    # Read TSV
    df = pd.read_csv(tsv_path, delimiter=delimiter)

    # Validate columns
    required_cols = [barcode_col, true_col]
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")

    # Create output DataFrame
    result_df = pd.DataFrame({
        'barcode': df[barcode_col],
        'annotation': df[true_col]
    })

    # Save to TSV
    output_path = os.path.join(out_dir, 'ground_truth.tsv')
    os.makedirs(out_dir, exist_ok=True)
    result_df.to_csv(output_path, sep='\t', index=False)
    print(f"Ground truth saved to {output_path}")

    # Print summary
    n_cells = len(result_df)
    n_unique_labels = len(result_df['annotation'].unique())
    print(f"Extracted {n_cells} cells with ground truth annotations.")
    print(f"Unique labels in {true_col}: {result_df['annotation'].unique()}")

# cluster for xclone and infercnv
# bcd/tumor_nontumor/utils.py   (or wherever tumor_classify lives)
import os
import pandas as pd
import numpy as np
from anndata import AnnData
from scipy.cluster.hierarchy import linkage, fcluster

def tumor_classify(adata, layer_name):
    """
    Cluster cells based on expression matrix using hierarchical clustering with Ward linkage
    and Euclidean distance, forcing 2 clusters. Label the cluster with the lowest aberration score
    (mean absolute deviation from median expression) as 'normal' and the other as 'tumor'.
    Parameters:
    -----------
    adata : anndata.AnnData
    AnnData object containing the expression matrix in adata.layers[layer_name].
    layer_name : str
    Name of the layer containing the expression matrix (shape: n_cells, n_genes).
    Returns:
    --------
    anndata.AnnData
    Updated AnnData object with 'tumor_pred' column in obs ('normal', 'tumor').
    """
    # Step 1: Extract the expression matrix
    if layer_name == "infercnv":
        expr_matrix = adata.X
    else:
        expr_matrix = adata.layers[layer_name] # Shape: (n_cells, n_genes)
    
    # Verify matrix shape
    n_cells, n_genes = expr_matrix.shape
    if expr_matrix.shape[0] != adata.n_obs:
        raise ValueError("Cell count mismatch between matrix and obs")

    # Step 2: Apply hierarchical clustering with Ward linkage and Euclidean distance
    Z = linkage(expr_matrix, method='ward', metric='euclidean')
    cluster_labels = fcluster(Z, t=2, criterion='maxclust') - 1 # Labels: 0 or 1
    n_clusters = 2

    # Step 3: Assign cluster labels based on aberration score
    # Compute aberration score as mean absolute deviation from median expression per gene
    if layer_name == "infercnv":
        # median_expr = np.median(expr_matrix, axis=0) # Median per gene
        aberration_scores = np.mean(np.abs(expr_matrix - 1), axis=1) # Shape: (n_cells,)
    else: #xclone
        aberration_scores = np.mean(np.abs(expr_matrix - 0), axis=1) # xclone: processed log ratio
    
    # Compute mean aberration scores per cluster
    mean_scores = [np.mean(aberration_scores[cluster_labels == i]) for i in [0, 1]]
    normal_cluster = np.argmin(mean_scores)
    tumor_pred = np.where(cluster_labels == normal_cluster, 'normal', 'tumor')

    # Step 4: Add predictions to adata.obs
    adata.obs['tumor_pred'] = tumor_pred
    # Print summary
    print(f"Processed {n_cells} cells and {n_genes} genes.")
    print(f"Mean aberration scores per cluster: {mean_scores}")
    print(f"Cluster labels: {['normal', 'tumor']}")
    for label in np.unique(tumor_pred):
        count = np.sum(tumor_pred == label)
    print(f"Number of cells in {label}: {count}")
    return adata


def xclone_predict_and_save(
    adata_path,
    layer_name,
    output_dir='output'
):
    """
    Read AnnData, perform tumor classification using hierarchical clustering, and save two TSV files:
    1. Ground truth (barcode, annotation) from true_col.
    2. Predictions (barcode, prediction) from tumor_classify.

    Parameters:
    -----------
    adata_path : str
        Path to .h5ad file containing the expression matrix and true annotations.
    layer_name : str
        Name of the layer containing the expression matrix (shape: n_cells, n_genes).
    output_dir : str
        Directory to save ground truth and predictions TSV files (default: 'output').

    Returns:
    --------
    None
        Saves two TSV files: ground_truth.tsv and xclone_predictions.tsv.
    """
    # Validate inputs
    if not os.path.exists(adata_path):
        raise ValueError(f"AnnData file not found at {adata_path}")

    # Read AnnData
    adata = ad.read_h5ad(adata_path)
    if layer_name not in adata.layers:
        raise ValueError(f"Layer '{layer_name}' not found in adata.layers")

    # Perform tumor classification
    adata = tumor_classify(adata, layer_name)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Save predictions TSV
    predictions_df = pd.DataFrame({
        'barcode': adata.obs_names,
        'prediction': adata.obs['tumor_pred']
    })
    predictions_path = os.path.join(output_dir, 'xclone_predictions.tsv')
    predictions_df.to_csv(predictions_path, sep='\t', index=False)
    print(f"Predictions saved to {predictions_path}")

# infercnv
def infercnv_predict_and_save(
    adata_path,
    output_dir='output'
):
    """
    Read AnnData, perform tumor classification using hierarchical clustering on adata.X,
    and save two TSV files: ground truth (barcode, annotation) if true_col is provided,
    and predictions (barcode, prediction).

    Parameters:
    -----------
    adata_path : str
        Path to .h5ad file containing the expression matrix in adata.X.
    output_dir : str
        Directory to save ground truth (if true_col provided) and predictions TSV files (default: 'output').
    true_col : str, optional
        adata.obs column with true annotations (default: None, skips ground truth TSV).

    Returns:
    --------
    None
        Saves TSV files: ground_truth.tsv (if true_col provided) and infercnv_predictions.tsv.
    """
    # Validate inputs
    if not os.path.exists(adata_path):
        raise ValueError(f"AnnData file not found at {adata_path}")

    # Read AnnData
    adata = ad.read_h5ad(adata_path)

    # Perform tumor classification
    adata = tumor_classify(adata, 
        layer_name="infercnv"
    )

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)


    # Save predictions TSV
    predictions_df = pd.DataFrame({
        'barcode': adata.obs['cell'],
        'prediction': adata.obs['tumor_pred']
    })
    predictions_path = os.path.join(output_dir, 'infercnv_predictions.tsv')
    predictions_df.to_csv(predictions_path, sep='\t', index=False)
    print(f"Predictions saved to {predictions_path}")

# numbat
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
import os

def numbat_process_predictions(
    tsv_path,
    out_dir='output',
    barcode_col='cell',
    p_cnv_col='p_cnv',
    delimiter='\t',
    n_clusters=2,
    random_state=42
):
    """
    Read Numbat TSV file, apply K-means clustering on p_cnv column to classify cells as 'normal' or 'tumor',
    and save predictions to a TSV file in out_dir.

    Parameters:
    -----------
    tsv_path : str
        Path to Numbat TSV file with columns including barcode_col and p_cnv_col.
    out_dir : str
        Directory to save the output TSV file (numbat_predictions.tsv).
    barcode_col : str
        Column name for cell barcodes (default: 'cell').
    p_cnv_col : str
        Column name for CNV probabilities (default: 'p_cnv').
    delimiter : str
        Delimiter for the input TSV file (default: '\t').
    n_clusters : int
        Number of clusters for K-means (default: 2 for normal/tumor).
    random_state : int
        Random seed for K-means (default: 42).

    Returns:
    --------
    None
        Saves a TSV file with columns: barcode, prediction ('normal' or 'tumor') to out_dir/numbat_predictions.tsv.
    """
    # Validate input file
    if not os.path.exists(tsv_path):
        raise ValueError(f"TSV file not found at {tsv_path}")

    # Read TSV
    df = pd.read_csv(tsv_path, delimiter=delimiter)

    # Validate columns
    required_cols = [barcode_col, p_cnv_col]
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")

    # Filter out rows with empty p_cnv
    initial_n_cells = len(df)
    df = df.dropna(subset=[p_cnv_col])
    n_cells = len(df)
    n_removed = initial_n_cells - n_cells
    if n_removed > 0:
        print(f"Removed {n_removed} cells with empty {p_cnv_col} values.")

    # Check if any cells remain
    if n_cells == 0:
        raise ValueError(f"No cells remain after removing empty {p_cnv_col} values.")

    # Extract p_cnv values
    p_cnv_values = df[p_cnv_col].values.reshape(-1, 1)

    # Apply K-means clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state)
    cluster_labels = kmeans.fit_predict(p_cnv_values)

    # Identify tumor cluster (higher mean p_cnv)
    mean_scores = np.zeros(n_clusters)
    for cluster in range(n_clusters):
        cluster_cells = cluster_labels == cluster
        mean_scores[cluster] = np.mean(p_cnv_values[cluster_cells])
    tumor_cluster = np.argmax(mean_scores)
    predictions = np.where(cluster_labels == tumor_cluster, 'tumor', 'normal')

    # Compute threshold (midpoint between cluster centers)
    cluster_centers = kmeans.cluster_centers_.flatten()
    low_center, high_center = sorted(cluster_centers)
    threshold = (low_center + high_center) / 2

    # Create output DataFrame
    result_df = pd.DataFrame({
        'barcode': df[barcode_col],
        'prediction': predictions
    })

    # Save to TSV
    output_path = os.path.join(out_dir, 'numbat_predictions.tsv')
    os.makedirs(out_dir, exist_ok=True)
    result_df.to_csv(output_path, sep='\t', index=False)
    print(f"Predictions saved to {output_path}")

    # Print summary
    n_tumor = np.sum(predictions == 'tumor')
    print(f"Processed {n_cells} cells after filtering.")
    print(f"Numbat p_cnv cluster centers: {cluster_centers}")
    print(f"Selected threshold: {threshold:.4f} (cells > threshold classified as tumor)")
    print(f"Number of tumor cells: {n_tumor}")
    print(f"Number of normal cells: {n_cells - n_tumor}")

# calicost
import pandas as pd
import os

def copykat_process_predictions(
    tsv_path,
    output_dir='output',
    delimiter='\t'
):
    """
    Read CopyKAT TSV file, rename columns to 'barcode' and 'prediction', convert 'diploid' to 'normal'
    and 'aneuploid' to 'tumor', and save to a new TSV file.

    Parameters:
    -----------
    tsv_path : str
        Path to CopyKAT TSV file with columns 'cell.names' and 'copykat.pred'.
    output_dir : str
        Dir to save the processed TSV file.
    delimiter : str
        Delimiter for the input TSV file (default: '\t').

    Returns:
    --------
    None
        Saves a TSV file with columns: barcode, prediction ('normal' or 'tumor').
    """
    # Validate input file
    if not os.path.exists(tsv_path):
        raise ValueError(f"TSV file not found at {tsv_path}")

    # Read TSV
    df = pd.read_csv(tsv_path, delimiter=delimiter)

    # Validate columns
    required_cols = ['cell.names', 'copykat.pred']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")

    # Validate prediction values
    valid_preds = {'diploid', 'aneuploid'}
    if not set(df['copykat.pred']).issubset(valid_preds):
        invalid_preds = set(df['copykat.pred']) - valid_preds
        raise ValueError(f"Invalid prediction values found: {invalid_preds}. Expected: {valid_preds}")

    # Create output DataFrame
    result_df = pd.DataFrame({
        'barcode': df['cell.names'],
        'prediction': df['copykat.pred'].replace({'diploid': 'normal', 'aneuploid': 'tumor'})
    })

    # Save to TSV
    predictions_path = os.path.join(output_dir, 'copykat_predictions.tsv')
    result_df.to_csv(predictions_path, sep='\t', index=False)
    print(f"Processed predictions saved to {output_dir}")

    # Print summary
    n_cells = len(df)
    n_tumor = sum(result_df['prediction'] == 'tumor')
    print(f"Processed {n_cells} cells.")
    print(f"Number of tumor cells: {n_tumor}")
    print(f"Number of normal cells: {n_cells - n_tumor}")

# calicost
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
import os

def predict_tumor_calicost(
    tsv_path,
    proportion_col='tumor_proportion',
    output_dir='output',
    delimiter='\t',
    n_clusters=2,
    random_state=42
):
    """
    Predict tumor cells from CalicoST output using K-means on tumor_proportion, excluding cells with empty tumor_proportion,
    and save to a TSV file with columns 'barcode' and 'prediction' ('tumor' or 'normal').

    Parameters:
    -----------
    tsv_path : str
        Path to CalicoST TSV file (columns: BARCODES, clone_label, tumor_proportion).
    proportion_col : str
        Column name for tumor proportions (default: 'tumor_proportion').
    output_dir : str
        Dir to save output TSV file.
    delimiter : str
        Delimiter for TSV file (default: '\t').
    n_clusters : int
        Number of clusters for K-means (default: 2 for tumor/non-tumor).
    random_state : int
        Random seed for K-means (default: 42).

    Returns:
    --------
    None
        Saves a TSV file with columns: barcode, prediction ('tumor' or 'normal') for cells with valid tumor_proportion.
    """
    # Validate output_fir
    if not output_dir:
        raise ValueError("output_dir must be provided to save predictions.")

    # Read CalicoST TSV
    df = pd.read_csv(tsv_path, delimiter=delimiter)

    # Validate columns
    required_cols = ['BARCODES', 'clone_label', proportion_col]
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")

    # Filter out rows with empty tumor_proportion
    initial_n_cells = len(df)
    df = df.dropna(subset=[proportion_col])
    n_cells = len(df)
    n_removed = initial_n_cells - n_cells
    if n_removed > 0:
        print(f"Removed {n_removed} cells with empty {proportion_col} values.")

    # Check if any cells remain
    if n_cells == 0:
        raise ValueError(f"No cells remain after removing empty {proportion_col} values.")

    # Extract tumor proportions
    tumor_proportions = df[proportion_col].values.reshape(-1, 1)

    # Apply K-means clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state)
    cluster_labels = kmeans.fit_predict(tumor_proportions)

    # Identify tumor cluster (higher mean tumor_proportion)
    mean_scores = np.zeros(n_clusters)
    for cluster in range(n_clusters):
        cluster_cells = cluster_labels == cluster
        mean_scores[cluster] = np.mean(tumor_proportions[cluster_cells])
    tumor_cluster = np.argmax(mean_scores)
    predictions = np.where(cluster_labels == tumor_cluster, 'tumor', 'normal')

    # Compute threshold (midpoint between cluster centers)
    cluster_centers = kmeans.cluster_centers_.flatten()
    low_center, high_center = sorted(cluster_centers)
    threshold = (low_center + high_center) / 2

    # Create output DataFrame
    result_df = pd.DataFrame({
        'barcode': df['BARCODES'],
        'prediction': predictions
    })

    # Save to TSV
    os.makedirs(os.path.dirname(output_dir) or '.', exist_ok=True)
    predictions_path = os.path.join(output_dir, 'calicost_predictions.tsv')
    result_df.to_csv(predictions_path, sep='\t', index=False)
    print(f"Predictions saved to {output_dir}")

    # Print summary
    n_tumor = np.sum(predictions == 'tumor')
    print(f"Processed {n_cells} cells after filtering.")
    print(f"CalicoST tumor_proportion cluster centers: {cluster_centers}")
    print(f"Selected threshold: {threshold:.4f} (cells > threshold classified as tumor)")
    print(f"Number of tumor cells: {n_tumor}")
    print(f"Number of non-tumor cells: {n_cells - n_tumor}")







# merge data
import pandas as pd
import os

def merge_predictions_to_tsv(
    tools,
    out_dir,
    true_col='annotation',
    true_label='normal',
    output_tsv_path=None
):
    """
    Read ground truth and prediction TSV files for specified tools from out_dir, intersect barcodes,
    merge into a DataFrame with string labels 'normal' and 'tumor', and save to a TSV file.

    Parameters:
    -----------
    tools : list of str
        List of tool names (e.g., ['xclone', 'infercnv', 'copykat', 'calicost']).
    out_dir : str
        Directory containing ground_truth.tsv and <tool>_predictions.tsv files.
    true_col : str
        Column in ground_truth.tsv with true labels (default: 'annotation').
    true_label : str
        Value in true_col for non-tumor cells (default: 'normal').
    output_tsv_path : str, optional
        Path to save the combined DataFrame as TSV (default: out_dir/combined_predictions.tsv).

    Returns:
    --------
    pandas.DataFrame
        DataFrame with columns: barcode, true_label, <tool>_pred, all using 'normal' or 'tumor'.
    """
    # Set default output path
    if output_tsv_path is None:
        output_tsv_path = os.path.join(out_dir, 'combined_predictions.tsv')

    # Validate inputs
    if not os.path.exists(out_dir):
        raise ValueError(f"Output directory not found: {out_dir}")

    # Read ground truth
    ground_truth_path = os.path.join(out_dir, 'ground_truth.tsv')
    if not os.path.exists(ground_truth_path):
        raise ValueError(f"Ground truth file not found: {ground_truth_path}")
    ground_truth_df = pd.read_csv(ground_truth_path, sep='\t')
    if not all(col in ground_truth_df.columns for col in ['barcode', true_col]):
        raise ValueError(f"Ground truth TSV must contain columns: ['barcode', '{true_col}']")

    # Initialize set of barcodes
    intersect_barcodes = set(ground_truth_df['barcode'])

    # Dictionary to store predictions
    predictions = {}
    for tool in tools:
        tsv_path = os.path.join(out_dir, f'{tool}_predictions.tsv')
        if not os.path.exists(tsv_path):
            raise ValueError(f"Prediction file not found for {tool}: {tsv_path}")
        df = pd.read_csv(tsv_path, sep='\t')
        if not all(col in df.columns for col in ['barcode', 'prediction']):
            raise ValueError(f"TSV for {tool} must contain columns: ['barcode', 'prediction']")
        # Validate prediction values
        valid_preds = {'normal', 'tumor', 'benign'}
        if not set(df['prediction']).issubset(valid_preds):
            invalid_preds = set(df['prediction']) - valid_preds
            raise ValueError(f"Invalid prediction values for {tool}: {invalid_preds}. Expected: {valid_preds}")
        intersect_barcodes &= set(df['barcode'])
        predictions[tool] = df.set_index('barcode')['prediction']

    # Check for intersecting barcodes
    intersect_barcodes = sorted(intersect_barcodes)
    if not intersect_barcodes:
        raise ValueError("No intersecting barcodes found across ground truth and prediction files.")

    # Subset ground truth to intersecting barcodes
    ground_truth_df = ground_truth_df[ground_truth_df['barcode'].isin(intersect_barcodes)].set_index('barcode').loc[intersect_barcodes]
    # Map ground truth to 'normal' or 'tumor'
    true_labels = ground_truth_df[true_col].apply(lambda x: 'normal' if x == true_label else 'tumor')

    # Create DataFrame
    result_df = pd.DataFrame({
        'barcode': intersect_barcodes,
        'true_label': true_labels.values
    })

    # Add predictions for each tool
    for tool in tools:
        pred_series = predictions[tool].loc[intersect_barcodes]
        # Map 'benign' to 'normal', keep 'normal' and 'tumor' as is
        result_df[f'{tool}_pred'] = pred_series.apply(lambda x: 'normal' if x in ['normal', 'benign'] else 'tumor').values

    # Save DataFrame to TSV
    os.makedirs(os.path.dirname(output_tsv_path) or '.', exist_ok=True)
    result_df.to_csv(output_tsv_path, sep='\t', index=False)
    print(f"Combined predictions saved to {output_tsv_path}")
    print(f"Processed {len(intersect_barcodes)} intersecting barcodes.")

    # Print summary
    n_true_non_tumor = sum(true_labels == 'normal')
    print(f"True non-tumor cells ({true_col} == {true_label}): {n_true_non_tumor}")
    for tool in tools:
        n_pred_tumor = sum(result_df[f'{tool}_pred'] == 'tumor')
        print(f"{tool} tumor cells: {n_pred_tumor}")

    return result_df


# evaluate
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score,
    f1_score, adjusted_rand_score, confusion_matrix
)


def evaluate_predictions(
    input_df,
    output_dir='output',
    true_label_col='true_label',
    pred_prefix='_pred'
):
    """
    Compute accuracy, precision, recall, F1 score, and ARI for each tool's predictions,
    plot ONE radar plot (all metrics), ONE confusion-matrix panel, and ONE histogram,
    and generate confusion matrix heatmaps.
    Handles string labels 'normal' and 'tumor'.

    Parameters
    ----------
    input_df : pandas.DataFrame or str
        DataFrame or path to TSV with columns: barcode, true_label, <method>_pred.
    output_dir : str
        Directory to save the three grouped PNG files.
    true_label_col : str
        Column with true labels (default: 'true_label').
    pred_prefix : str
        Suffix of prediction columns (default: '_pred').

    Returns
    -------
    dict
        {'accuracy':{...}, 'precision':{...}, 'recall':{...}, 'f1':{...}, 'ari':{...}}
    """
    # ------------------------------------------------------------------ #
    # 1. Load data
    # ------------------------------------------------------------------ #
    if isinstance(input_df, str):
        df = pd.read_csv(input_df, sep='\t')
    else:
        df = input_df.copy()

    # ------------------------------------------------------------------ #
    # 2. Validation
    # ------------------------------------------------------------------ #
    if true_label_col not in df.columns:
        raise ValueError(f"Column '{true_label_col}' not found.")
    pred_cols = [c for c in df.columns if c.endswith(pred_prefix)]
    if not pred_cols:
        raise ValueError(f"No columns ending with '{pred_prefix}'.")
    valid = {'normal', 'tumor'}
    if not set(df[true_label_col]).issubset(valid):
        raise ValueError(f"Invalid true labels: {set(df[true_label_col]) - valid}")
    for c in pred_cols:
        if not set(df[c]).issubset(valid):
            raise ValueError(f"Invalid labels in {c}")

    # ------------------------------------------------------------------ #
    # 3. Initialise metric containers (exactly the original set)
    # ------------------------------------------------------------------ #
    methods = [c.replace(pred_prefix, '') for c in pred_cols]
    metrics = {
        'accuracy': {}, 'precision': {}, 'recall': {}, 'f1': {}, 'ari': {}
    }

    # ------------------------------------------------------------------ #
    # 4. Compute metrics
    # ------------------------------------------------------------------ #
    y_true = df[true_label_col] == 'tumor'                 # boolean for tumor
    for pred_col in pred_cols:
        method = pred_col.replace(pred_prefix, '')
        y_pred = df[pred_col] == 'tumor'

        metrics['accuracy'][method] = accuracy_score(y_true, y_pred)
        metrics['precision'][method] = precision_score(
            y_true, y_pred, zero_division=0
        )
        metrics['recall'][method] = recall_score(
            y_true, y_pred, zero_division=0
        )
        metrics['f1'][method] = f1_score(y_true, y_pred, zero_division=0)
        metrics['ari'][method] = adjusted_rand_score(
            df[true_label_col], df[pred_col]
        )

    # ------------------------------------------------------------------ #
    # 5. Create output directory
    # ------------------------------------------------------------------ #
    os.makedirs(output_dir, exist_ok=True)

    # ------------------------------------------------------------------ #
    # 6. 1. CONFUSION-MATRIX PANEL (all tools in ONE figure)
    # ------------------------------------------------------------------ #
    n_tools = len(methods)
    fig_cm, axes_cm = plt.subplots(
        1, n_tools, figsize=(4 * n_tools, 4)
    )
    if n_tools == 1:
        axes_cm = [axes_cm]

    for ax, method, pred_col in zip(axes_cm, methods, pred_cols):
        cm = confusion_matrix(
            df[true_label_col], df[pred_col], labels=['normal', 'tumor']
        )
        sns.heatmap(
            cm, annot=True, fmt='d', cmap='Blues',
            xticklabels=['normal', 'tumor'],
            yticklabels=['normal', 'tumor'],
            cbar=False, ax=ax
        )
        ax.set_title(f'{method}') # \nARI={metrics["ari"][method]:.3f}
        ax.set_xlabel('Predicted')
        ax.set_ylabel('True')
    plt.tight_layout()
    cm_path = os.path.join(output_dir, 'confusion_matrices.png')
    fig_cm.savefig(cm_path, dpi=300, bbox_inches='tight')
    plt.close(fig_cm)
    print(f"Confusion-matrix panel → {cm_path}")

    # ------------------------------------------------------------------ #
    # 6. 2. RADAR PANEL (single polar chart, all 5 metrics)
    # ------------------------------------------------------------------ #
    metric_names = ['accuracy', 'precision', 'recall', 'f1', 'ari']
    angles = np.linspace(0, 2*np.pi, len(metric_names), endpoint=False).tolist()
    angles += angles[:1]                     # close circle

    fig_r, ax_r = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    for method in methods:
        values = [metrics[m][method] for m in metric_names]
        values += values[:1]
        ax_r.plot(angles, values, 'o-', linewidth=2, label=method)
        ax_r.fill(angles, values, alpha=0.1)
    ax_r.set_xticks(angles[:-1])
    ax_r.set_xticklabels([m.title() for m in metric_names])
    ax_r.set_ylim(0, 1)
    ax_r.set_title('Metric Comparison (Radar)', size=16, pad=20)
    ax_r.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))

    radar_path = os.path.join(output_dir, 'metrics_radar.png')
    fig_r.savefig(radar_path, dpi=300, bbox_inches='tight')
    plt.close(fig_r)
    print(f"Radar panel → {radar_path}")

    # ------------------------------------------------------------------ #
    # 6. 3. HISTOGRAM PANEL (predicted label counts)
    # ------------------------------------------------------------------ #
    fig_h, ax_h = plt.subplots(figsize=(10, 5))
    bar_width = 0.35
    pos = np.arange(len(methods))

    tumor_counts = [(df[f'{m}{pred_prefix}'] == 'tumor').sum() for m in methods]
    normal_counts = [len(df) - t for t in tumor_counts]

    ax_h.bar(pos - bar_width/2, normal_counts, bar_width,
             label='normal', color='#1f77b4')
    ax_h.bar(pos + bar_width/2, tumor_counts, bar_width,
             label='tumor', color='#ff7f0e')

    ax_h.set_xticks(pos)
    ax_h.set_xticklabels(methods)
    ax_h.set_ylabel('Cell count')
    ax_h.set_title('Predicted Label Distribution')
    ax_h.legend()

    # percentages on top of bars
    total = len(df)
    for i, (n, t) in enumerate(zip(normal_counts, tumor_counts)):
        ax_h.text(i - bar_width/2, n + total*0.01, f'{n/total:.1%}',
                 ha='center', va='bottom', fontsize=9)
        ax_h.text(i + bar_width/2, t + total*0.01, f'{t/total:.1%}',
                 ha='center', va='bottom', fontsize=9)

    hist_path = os.path.join(output_dir, 'prediction_histogram.png')
    fig_h.savefig(hist_path, dpi=300, bbox_inches='tight')
    plt.close(fig_h)
    print(f"Histogram panel → {hist_path}")

    # ------------------------------------------------------------------ #
    # 6. 4. METRIC BAR-PLOT PANEL (one subplot per metric)
    # ------------------------------------------------------------------ #
    metric_names = ['accuracy', 'precision', 'recall', 'f1', 'ari']
    n_metrics    = len(metric_names)

    # Create a grid: 1 row × n_metrics columns (or 2 rows if many metrics)
    fig_b, axes_b = plt.subplots(
        1, n_metrics, figsize=(4 * n_metrics, 5), sharey=True
    )
    if n_metrics == 1:
        axes_b = [axes_b]

    bar_width = 0.6
    pos       = np.arange(len(methods))

    for ax, metric in zip(axes_b, metric_names):
        values = [metrics[metric][m] for m in methods]

        bars = ax.bar(pos, values, bar_width, color='#1f77b4', edgecolor='black')
        ax.set_title(metric.title(), fontsize=12, fontweight='bold')
        ax.set_xticks(pos)
        ax.set_xticklabels(methods, rotation=45, ha='right')
        ax.set_ylim(0, 1)

        # print value on top of each bar
        for bar, val in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.02,
                f'{val:.3f}',
                ha='center', va='bottom', fontsize=9
            )

    # common y-label
    fig_b.supylabel('Score', fontsize=12)
    fig_b.suptitle('Metric Comparison (Bar Plot)', fontsize=14, y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    bar_path = os.path.join(output_dir, 'metrics_barplot.png')
    fig_b.savefig(bar_path, dpi=300, bbox_inches='tight')
    plt.close(fig_b)
    print(f"Metric bar-plot panel → {bar_path}")

    # ------------------------------------------------------------------ #
    # 7. SAVE METRICS TO TSV
    # ------------------------------------------------------------------ #
    # Build a DataFrame: rows = methods, columns = metrics
    metrics_df = pd.DataFrame(metrics)
    metrics_df = metrics_df.reindex(columns=metric_names)  # ensure order
    metrics_df.index.name = 'method'
    tsv_path = os.path.join(output_dir, 'metrics_summary.tsv')
    metrics_df.to_csv(tsv_path, sep='\t')
    print(f"Metrics saved to TSV → {tsv_path}")

    # ------------------------------------------------------------------ #
    # 8. Console summary (original style)
    # ------------------------------------------------------------------ #
    print("\n=== Evaluation Summary ===")
    print(f"Cells evaluated : {len(df)}")
    print(f"True tumor cells: {sum(df[true_label_col]=='tumor')} "
          f"({sum(df[true_label_col]=='tumor')/len(df):.1%})")
    for m in methods:
        pred_tumor = sum(df[f'{m}{pred_prefix}'] == 'tumor')
        print(f"\n{m}:")
        print(f"  pred tumor  : {pred_tumor} ({pred_tumor/len(df):.1%})")
        print(f"  Accuracy    : {metrics['accuracy'][m]:.4f}")
        print(f"  Precision   : {metrics['precision'][m]:.4f}")
        print(f"  Recall      : {metrics['recall'][m]:.4f}")
        print(f"  F1          : {metrics['f1'][m]:.4f}")
        print(f"  ARI         : {metrics['ari'][m]:.4f}")

    return metrics