# metric.py - functions for calculating metrics for binary classification.



import gc
import numpy as np
import os
import pandas as pd
import scipy as sp
from logging import info
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score,
    f1_score, adjusted_rand_score
)
from ..utils.base import assert_e



def run_metric(
    tool_list,
    tool_fn_list,
    truth_fn,
    out_fn,
    verbose = True
):
    """Main function of calculating metrics.
    
    Parameters
    ----------
    tool_list : list of Tool
        A list of tool-specific :class:`~.tool.Tool` objects.
    tool_fn_list : list of str
        A list of tool-specific files storing predicted labels.
    truth_fn : str
        The file storing the ground truth labels.
    out_fn : str
        The output metrics file.
    verbose : bool, default True
        Whether to show detailed logging information.
        
    Returns
    -------
    dict
        Results.
    """
    # check args.
    assert len(tool_list) == len(tool_fn_list)
    for fn in tool_fn_list:
        assert_e(fn)
    assert_e(truth_fn)

    
    # calculate metrics for each tool.
    if verbose:
        info("calculate metrics for each tool ...")

    accuracy_list = []
    precision_list = []
    recall_list = []
    f1_list = []
    ari_list = []
    
    df = pd.read_csv(truth_fn, sep = '\t')
    truth_labels = df['annotation'].to_numpy()
    for i, (tool, tool_fn) in enumerate(zip(tool_list, tool_fn_list)):
        tid = tool.tid
        if verbose:
            info("process %s ..." % tid)
        
        df = pd.read_csv(tool_fn, sep = '\t')
        tool_labels = df['prediction'].to_numpy()
        res = calc_binary_metrics(
            truth = truth_labels,
            pred = tool_labels,
            pos_label = 'tumor'
        )
        accuracy_list.append(res['accuracy'])
        precision_list.append(res['precision'])
        recall_list.append(res['recall'])
        f1_list.append(res['F1'])
        ari_list.append(res['ARI'])


    # save files.
    if verbose:
        info("save metric result files ...")

    df_metric = pd.DataFrame(
        data = dict(
            tool = [tool.tid for tool in tool_list],
            accuracy = accuracy_list,
            precision = precision_list,
            recall = recall_list,
            F1 = f1_list,
            ARI = ari_list
        ))
    df_metric = pd.melt(
        df_metric, 
        id_vars = ['tool'],
        var_name = 'metric', value_name = 'value'
    )
    df_metric.to_csv(out_fn, sep = "\t", index = False)
    
    
    if verbose:
        info("Summary of metrics:")
        info(str(df_metric))


    res = dict(
        # out_fn : str
        #   File storing DataFrame that contains multiple metrics.
        #   It has six columns "tool", "accuracy", "precision", "recall",
        #   "F1", and "ARI".
        out_fn = out_fn
    )
    return(res)



def calc_binary_metrics(truth, pred, pos_label):
    """
    truth : array
        Ground truth labels.
    pred : array
        Predicted labels.
    pos_label
        Positive label.
    """
    y_pred = (np.array(pred) == pos_label) + 0
    y_true = (np.array(truth) == pos_label) + 0

    res = dict(
        accuracy = accuracy_score(y_true, y_pred),
        precision = precision_score(
            y_true, y_pred, zero_division = 0
        ),
        recall = recall_score(
            y_true, y_pred, zero_division = 0
        ),
        F1 = f1_score(y_true, y_pred, zero_division = 0),
        ARI = adjusted_rand_score(y_true, y_pred)
    )
    return(res)

