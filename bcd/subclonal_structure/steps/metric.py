# metric.py



import gc
import numpy as np
import os
import pandas as pd
import scipy as sp
from logging import info
from sklearn.metrics import adjusted_rand_score
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

    ari_list = []

    df = pd.read_csv(truth_fn, sep = '\t')
    truth_labels = df['annotation'].to_numpy()
    for i, (tool, tool_fn) in enumerate(zip(tool_list, tool_fn_list)):
        tid = tool.tid
        if verbose:
            info("process %s ..." % tid)
        
        df = pd.read_csv(tool_fn, sep = '\t')
        tool_labels = df['prediction'].to_numpy()
        res = calc_metrics(
            truth = truth_labels,
            pred = tool_labels
        )
        ari_list.append(res['ARI'])


    # save files.
    if verbose:
        info("save metric result files ...")

    df_metric = pd.DataFrame(
        data = dict(
            tool = [tool.tid for tool in tool_list],
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
        #   It has two columns "tool", and "ARI".
        out_fn = out_fn
    )
    return(res)



def calc_metrics(truth, pred):
    """
    truth : array
        Ground truth labels.
    pred : array
        Predicted labels.
    """
    y_pred = np.array(pred)
    y_true = np.array(truth)

    res = dict(
        ARI = adjusted_rand_score(y_true, y_pred)
    )
    return(res)
