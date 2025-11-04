# Tumor non_tumor classification benchmark main
# Author: Jiamu James Qiao
# Date: 2025-10

import os
from .utils import (
    extract_ground_truth,
    xclone_predict_and_save,
    infercnv_predict_and_save,
    copykat_process_predictions,
    predict_tumor_calicost,
    numbat_process_predictions,
    merge_predictions_to_tsv,
    evaluate_predictions
)


"""
Example usage:

tool_configs = [
    ('xclone', 'path/to/your/file.h5ad', {'layer_name': 'layer_name'}), # default: 'WMA_smoothed_log_ratio_ab_dynamic'
    ('infercnv', 'path/to/your/file.h5ad'),
    ('copykat', 'path/to/copykat_predictions.txt', {'delimiter': '\t'}),
    ('calicost', 'path/to/clone_labels.tsv', {'proportion_col': 'tumor_proportion', 'delimiter': '\t'})
]
tool_configs = [
    ('xclone', '$ADATA_PATH', {'layer_name': '$LAYER_NAME'}),
    ('infercnv', '$ADATA_PATH'),
    ('copykat', '$COPYKAT_TSV', {'delimiter': '\t'}),
    ('calicost', '$CALICOST_TSV', {'proportion_col': 'tumor_proportion', 'delimiter': '\t'}),
    ('numbat', '$NUMBAT_TSV', {'barcode_col': 'cell', 'p_cnv_col': 'p_cnv', 'delimiter': '\t'})
]
metrics = run_tumor_prediction_pipeline(
    tool_configs=tool_configs,
    out_dir='output',
    true_col='annotation',
    true_label='normal'
)
"""

def run_tumor_prediction_pipeline(
    tool_configs,
    ground_truth_tsv,
    out_dir='output',
    true_col='spot_anno',
    true_label='normal',
    barcode_col='barcode',
    gt_delimiter='\t'
):
    """
    Run the tumor prediction and evaluation pipeline for specified tools, merging and evaluating predictions.

    Parameters:
    -----------
    tool_configs : list of tuples
        List of (tool, input_path, **kwargs) where:
        - tool: str, tool name ('xclone', 'infercnv', 'copykat', 'calicost', 'numbat').
        - input_path: str, path to input file (.h5ad for xclone/infercnv, TSV for copykat/calicost/numbat).
        - **kwargs: tool-specific parameters (e.g., layer_name for xclone, proportion_col for calicost).
    ground_truth_tsv : str
        Path to TSV file containing ground truth annotations.
    out_dir : str
        Directory to save output TSVs and plots.
    true_col : str
        Column in ground truth TSV with true labels (default: 'annotation').
    true_label : str
        Value in true_col for non-tumor cells (default: 'normal').
    barcode_col : str
        Column in ground truth TSV with cell barcodes (default: 'barcode').
    gt_delimiter : str
        Delimiter for ground truth TSV (default: '\t').

    Returns:
    --------
    dict
        Dictionary with metrics (accuracy, precision, recall, f1, ari) for each tool.
    """
    # Validate inputs
    valid_tools = {'xclone', 'infercnv', 'copykat', 'calicost', 'numbat'}
    tools = []
    for config in tool_configs:
        if len(config) < 2:
            raise ValueError("Each tool_config must include tool name and input_path")
        tool = config[0].lower()
        if tool not in valid_tools:
            raise ValueError(f"Invalid tool: {tool}. Expected: {valid_tools}")
        input_path = config[1]
        if not os.path.exists(input_path):
            raise ValueError(f"Input file not found: {input_path}")
        tools.append(tool)

    # Validate ground truth TSV
    if not os.path.exists(ground_truth_tsv):
        raise ValueError(f"Ground truth TSV not found at {ground_truth_tsv}")

    # Create output directory
    os.makedirs(out_dir, exist_ok=True)

    # Step 1: Extract ground truth
    extract_ground_truth(
        tsv_path=ground_truth_tsv,
        out_dir=out_dir,
        barcode_col=barcode_col,
        true_col=true_col,
        delimiter=gt_delimiter
    )

    # Step 2: Generate predictions for each tool
    for config in tool_configs:
        tool, input_path, kwargs = config[0].lower(), config[1], config[2] if len(config) > 2 else {}
        
        if tool == 'xclone':
            layer_name = kwargs.get('layer_name', 'layer_name')
            xclone_predict_and_save(
                adata_path=input_path,
                layer_name=layer_name,
                output_dir=out_dir
            )
        elif tool == 'infercnv':
            infercnv_predict_and_save(
                rds_path=input_path,
                output_dir=out_dir
            )
        elif tool == 'copykat':
            copykat_process_predictions(
                tsv_path=input_path,
                output_dir=out_dir,
                delimiter=kwargs.get('delimiter', '\t')
            )
        elif tool == 'calicost':
            predict_tumor_calicost(
                tsv_path=input_path,
                proportion_col=kwargs.get('proportion_col', 'tumor_proportion'),
                output_dir=out_dir,
                delimiter=kwargs.get('delimiter', '\t')
            )
        elif tool == 'numbat':
            numbat_process_predictions(
                tsv_path=input_path,
                out_dir=out_dir,
                barcode_col=kwargs.get('barcode_col', 'cell'),
                p_cnv_col=kwargs.get('p_cnv_col', 'p_cnv'),
                delimiter=kwargs.get('delimiter', '\t'),
                n_clusters=kwargs.get('n_clusters', 2),
                random_state=kwargs.get('random_state', 42)
            )

    # Step 3: Merge predictions
    df = merge_predictions_to_tsv(
        tools=tools,
        out_dir=out_dir,
        true_col="annotation",
        true_label=true_label,
        output_tsv_path=os.path.join(out_dir, 'combined_predictions.tsv')
    )

    # Step 4: Evaluate predictions
    metrics = evaluate_predictions(
        input_df=df,
        output_dir=out_dir,
        true_label_col='true_label',
        pred_prefix='_pred'
    )

    print("Pipeline completed. Metrics:")
    print(metrics)
    return metrics

def merge_predict(
    tool_configs,
    ground_truth_tsv,
    out_dir='output',
    true_col='spot_anno',
    true_label='normal',
    barcode_col='barcode',
    gt_delimiter='\t'
):
    """
    Run the tumor prediction and evaluation pipeline for specified tools, merging and evaluating predictions.

    Parameters:
    -----------
    tool_configs : list of tuples
        List of (tool, input_path, **kwargs) where:
        - tool: str, tool name ('xclone', 'infercnv', 'copykat', 'calicost', 'numbat').
        - input_path: str, path to input file (.h5ad for xclone/infercnv, TSV for copykat/calicost/numbat).
        - **kwargs: tool-specific parameters (e.g., layer_name for xclone, proportion_col for calicost).
    ground_truth_tsv : str
        Path to TSV file containing ground truth annotations.
    out_dir : str
        Directory to save output TSVs and plots.
    true_col : str
        Column in ground truth TSV with true labels (default: 'annotation').
    true_label : str
        Value in true_col for non-tumor cells (default: 'normal').
    barcode_col : str
        Column in ground truth TSV with cell barcodes (default: 'barcode').
    gt_delimiter : str
        Delimiter for ground truth TSV (default: '\t').

    Returns:
    --------
    dict
        Dictionary with metrics (accuracy, precision, recall, f1, ari) for each tool.
    """
    # Validate inputs
    valid_tools = {'xclone', 'infercnv', 'copykat', 'calicost', 'numbat'}
    tools = []
    for config in tool_configs:
        if len(config) < 2:
            raise ValueError("Each tool_config must include tool name and input_path")
        tool = config[0].lower()
        if tool not in valid_tools:
            raise ValueError(f"Invalid tool: {tool}. Expected: {valid_tools}")
        input_path = config[1]
        if not os.path.exists(input_path):
            raise ValueError(f"Input file not found: {input_path}")
        tools.append(tool)

    # Validate ground truth TSV
    if not os.path.exists(ground_truth_tsv):
        raise ValueError(f"Ground truth TSV not found at {ground_truth_tsv}")

    # Create output directory
    os.makedirs(out_dir, exist_ok=True)

    '''
    # Step 1: Extract ground truth
    extract_ground_truth(
        tsv_path=ground_truth_tsv,
        out_dir=out_dir,
        barcode_col=barcode_col,
        true_col=true_col,
        delimiter=gt_delimiter
    )

    # Step 2: Generate predictions for each tool
    for config in tool_configs:
        tool, input_path, kwargs = config[0].lower(), config[1], config[2] if len(config) > 2 else {}
        
        if tool == 'xclone':
            layer_name = kwargs.get('layer_name', 'layer_name')
            xclone_predict_and_save(
                adata_path=input_path,
                layer_name=layer_name,
                output_dir=out_dir
            )
        elif tool == 'infercnv':
            infercnv_predict_and_save(
                rds_path=input_path,
                output_dir=out_dir
            )
        elif tool == 'copykat':
            copykat_process_predictions(
                tsv_path=input_path,
                output_dir=out_dir,
                delimiter=kwargs.get('delimiter', '\t')
            )
        elif tool == 'calicost':
            predict_tumor_calicost(
                tsv_path=input_path,
                proportion_col=kwargs.get('proportion_col', 'tumor_proportion'),
                output_dir=out_dir,
                delimiter=kwargs.get('delimiter', '\t')
            )
        elif tool == 'numbat':
            numbat_process_predictions(
                tsv_path=input_path,
                out_dir=out_dir,
                barcode_col=kwargs.get('barcode_col', 'cell'),
                p_cnv_col=kwargs.get('p_cnv_col', 'p_cnv'),
                delimiter=kwargs.get('delimiter', '\t'),
                n_clusters=kwargs.get('n_clusters', 2),
                random_state=kwargs.get('random_state', 42)
            )
    '''
    
    # Step 3: Merge predictions
    df = merge_predictions_to_tsv(
        tools=tools,
        out_dir=out_dir,
        true_col="annotation",
        true_label=true_label,
        output_tsv_path=os.path.join(out_dir, 'combined_predictions.tsv')
    )

    # Step 4: Evaluate predictions
    metrics = evaluate_predictions(
        input_df=df,
        output_dir=out_dir,
        true_label_col='true_label',
        pred_prefix='_pred'
    )

    print("Pipeline completed. Metrics:")
    print(metrics)
    return metrics