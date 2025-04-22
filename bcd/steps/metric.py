# metric.py - functions for calculating ROC and PRC.



import gc
import numpy as np
import os
import pandas as pd
from logging import info
from ..utils.base import assert_e
from ..utils.io import load_h5ad



def run_metric(
    args_list,
    out_dir,
    out_prefix,
    tool_fn_list,
    truth_fn,
    cna_type,
    max_n_cutoff = 1000,
    verbose = True
):
    """Main function of calculating metrics.
    
    Parameters
    ----------
    args_list : list of ToolArgs
        A list of tool-specific :class:`~.args.ToolArgs` objects.
    out_dir : str
        The output folder.
    out_prefix : str
        The prefix to output files.
    tool_fn_list : list of str
        A list of tool-specific files storing cell x gene expression or 
        probability matrix.
    truth_fn : str
        File storing the ground truth cell x gene expression or probability
        matrix.
    cna_type : str
        CNA types, one of {"gain", "loss", "loh"}.
    max_n_cutoff : int or None, default 1000
        Maximum number of cutoff values for calculating metrics.
        If None, use all unique values in tool matrix.
    verbose : bool, default True
        Whether to show detailed logging information.
        
    Returns
    -------
    dict
        Results.
    """
    # check args.
    assert len(args_list) == len(tool_fn_list)
    for fn in tool_fn_list:
        assert_e(fn)
    os.makedirs(out_dir, exist_ok = True)
    assert_e(truth_fn)
    assert cna_type in ("gain", "loss", "loh")

    
    # calculate metrics for each tool.
    if verbose:
        info("%s: calculate metrics for each tool ..." % cna_type)
    
    truth = load_h5ad(truth_fn)
    truth_mtx = truth.X
    
    df_metric = None
    auroc_list = []
    auprc_list = []
    for i, (args, tool_fn) in enumerate(zip(args_list, tool_fn_list)):
        tid = args.tid
        if verbose:
            info("process %s ..." % tid)
            
        adata = load_h5ad(tool_fn)
        tool_mtx = adata.X
        
        if tid.lower() == "infercnv" and cna_type == "loss":
            tool_mtx = tool_mtx * -1
        
        res = calc_metric(
            tid = tid,
            tool_mtx = tool_mtx,
            truth_mtx = truth_mtx,
            cna_type = cna_type,
            max_n_cutoff = max_n_cutoff,
            verbose = verbose
        )
        
        df = res["df"]
        df["tool"] = tid
        if i == 0:
            df_metric = df.copy()
        else:
            df_metric = pd.concat([df_metric, df], ignore_index = True)
        auroc_list.append(res["AUROC"])
        auprc_list.append(res["AUPRC"])
        
        del adata
        del df
        del res
        gc.collect()
        
        
    # save files.
    if verbose:
        info("save metric result files ...")
    
    roc_fn = os.path.join(out_dir, "%s.%s.roc.tsv" % \
                             (out_prefix, cna_type))
    df_roc = df_metric[["tool", "cutoff", "tpr", "fpr"]]
    df_roc.to_csv(roc_fn, sep = "\t", index = False)
    
    
    prc_fn = os.path.join(out_dir, "%s.%s.prc.tsv" % \
                             (out_prefix, cna_type))
    df_prc = df_metric[["tool", "cutoff", "precision"]].copy()
    df_prc["recall"] = df_metric["tpr"]
    df_prc.to_csv(prc_fn, sep = "\t", index = False)
    
    
    df_auroc = pd.DataFrame(
        data = dict(
            tool = [args.tid for args in args_list],
            AUROC = auroc_list
        ))
    auroc_fn = os.path.join(out_dir, "%s.%s.auroc.tsv" % \
                             (out_prefix, cna_type))
    df_auroc.to_csv(auroc_fn, sep = "\t", index = False)
    if verbose:
        info("%s: AUROC:\n\n%s\n" % (cna_type, str(df_auroc)))
    
    
    df_auprc = pd.DataFrame(
        data = dict(
            tool = [args.tid for args in args_list],
            AUPRC = auprc_list
        ))
    auprc_fn = os.path.join(out_dir, "%s.%s.auprc.tsv" % \
                             (out_prefix, cna_type))
    df_auprc.to_csv(auprc_fn, sep = "\t", index = False)
    if verbose:
        info("%s: AUPRC:\n\n%s\n" % (cna_type, str(df_auprc)))
        
        
    res = dict(
        # roc_fn : str
        #   File storing DataFrame that contains ROC at various cutoffs.
        #   It has four columns "tool", "cutoff", tpr", "fpr".
        roc_fn = roc_fn,
        
        # prc_fn : str
        #   File storing DataFrame that contains PRC at various cutoffs.
        #   It has four columns "tool", "cutoff", "precision", "recall".
        prc_fn = prc_fn,
        
        # auroc_fn : str
        #   File storing DataFrame that contains AUROC of each tool.
        #   It has two columns "tool" and "AUROC".
        auroc_fn = auroc_fn,
        
        # auprc_fn : str
        #   File storing DataFrame that contains AUPRC of each tool.
        #   It has two columns "tool" and "AUPRC".
        auprc_fn = auprc_fn
    )
    return(res)



def calc_metric(
    tid,
    tool_mtx,
    truth_mtx,
    cna_type,
    max_n_cutoff = 1000,
    verbose = True
):
    """Calculate ROC and PRC.
    
    Parameters
    ----------
    tid : str
        Tool ID.
    tool_mtx : numpy.matrix-like or scipy.sparse matrix
        The tool's cell x gene expression or probability matrix.
    truth_mtx : numpy.matrix-like or scipy.sparse matrix
        The ground truth cell x gene expression or probability matrix.
    cna_type : str
        CNA types, one of {"gain", "loss", "loh"}.
    max_n_cutoff : int or None, default 1000
        Maximum number of cutoff values for `metric`.
        If None, use all unique values in `tool_mtx`.
    verbose : bool, default True
        Whether to show detailed logging information.
        
    Returns
    -------
    dict
        Results.
    """
    # check args.
    assert np.all(np.array(tool_mtx.shape) == np.array(truth_mtx.shape))
    assert cna_type in ("gain", "loss", "loh")
    
    if verbose:
        info("%s-%s: tool matrix %s and truth matrix %s." % \
             (cna_type, tid, str(tool_mtx.shape), str(truth_mtx.shape)))
    
    tool_vec = tool_mtx.flatten()
    cutoff = np.unique(tool_vec)

    eps = (np.max(cutoff) - np.min(cutoff)) / len(cutoff)
    if eps <= 0:
        eps = 1e-6
    max_cutoff = np.max(cutoff) + eps
    min_cutoff = np.min(cutoff) - eps

    if max_n_cutoff is not None and len(cutoff) > max_n_cutoff - 2:
        cutoff = np.random.choice(cutoff, size = max_n_cutoff - 2, replace = False)
    cutoff = np.sort(np.unique(np.concatenate([cutoff, [min_cutoff, max_cutoff]])))
        
    if verbose:
        info("%s-%s: %d cutoffs." % (cna_type, tid, len(cutoff)))
    
    
    res = binary_metric(
        scores = tool_vec,
        labels = truth_mtx.toarray().flatten(),
        cutoff = cutoff,
        cut_direction = ">=",
        add_cut_inf = True
    )
        
    if verbose:
        info("%s-%s: AUROC = %f; AUPRC = %f;" % \
             (cna_type, tid, res["AUROC"], res["AUPRC"]))
    return(res)

        

# ref: the `binaryROC()` function in R package "cardelino".
# https://github.com/single-cell-genetics/cardelino/blob/3d2d9779877742dae8f3bc8644b5a0bee9973bd4/R/assessment.R#L192C1-L262C2
def binary_metric(
    scores,
    labels,
    cutoff = None,
    cut_direction = '>=',
    add_cut_inf = True
):
    """ROC and PRC curve for binary label prediction
    
    Parameters
    ----------
    scores : array-like of numeric
        Prediction score for each sample.
    labels : array-like of {0, 1}
        True labels for each sample.
    cutoff : array-like of numeric or None
        A vector of cutoffs; if None use all unique scores.
    cut_direction : {'>=', '>', '<=', '<'}
        A string to compare with cutoff.
    add_cut_inf : bool, default True
        If True, manually add a cutoff of `np.inf`.
        
    Returns
    -------
    dict
        Results.
    """
    # tpr = TP / (TP + FN)
    # fpr = FP / (FP + TN)
    # precision = TP / (TP + FP)
    # recall = tpr
    
    if cutoff is None:
        cutoff = np.sort(np.unique(scores))
    cutoff = np.sort(np.unique(np.concatenate([cutoff, [0, 1]])))

    n = len(cutoff)
    tpr = np.zeros(n)
    fpr = np.zeros(n)
    precision = np.zeros(n)
    
    idx = None
    for i in range(n):
        if cut_direction == "<=":
            idx = scores <= cutoff[i]
        elif cut_direction == "<":
            idx = scores < cutoff[i]
        elif cut_direction == ">":
            idx = scores > cutoff[i]
        else:
            idx = scores >= cutoff[i]
            
        if np.sum(idx) == 0:    # i.e., TP = FP = 0
            fpr[i] = 0.0
            tpr[i] = 0.0
            precision[i] = 1.0
        else:
            fpr[i] = np.sum(labels[idx] == 0) / np.sum(labels == 0)
            tpr[i] = np.sum(labels[idx] == 1) / np.sum(labels == 1)
            precision[i] = np.mean(labels[idx] == 1)
    #precision[logical_and(tpr == 0, np.isnan(precision))] = 1.0


    # Note, 
    # For ROC and PRC plot, it is important to set the leftmost and rightmost
    # of the x axis, which is `fpr` (in ROC plot) and `tpr` (in PRC plot),
    # respectively.
    # - `fpr` and `tpr` decreases when cutoff increases.
    # - The `fpr` and `tpr` can always reach (1.0, 1.0) by setting an extreamly
    #   small cutoff;
    # - However, it cannot guarantee they can always reach (0.0, 0.0) unless
    #   manually set.
    if add_cut_inf:
        cutoff = np.concatenate([cutoff, [np.inf]])
        fpr = np.concatenate([fpr, [0.0]])
        tpr = np.concatenate([tpr, [0.0]])
        precision = np.concatenate([precision, [1.0]])
        n += 1

    auroc = auprc = 0.0
    for i in range(n - 1):
        auroc += (fpr[i] - fpr[i+1]) * (tpr[i] + tpr[i+1]) * 0.5
        auprc += (tpr[i] - tpr[i+1]) * (precision[i] + precision[i+1]) * 0.5
    auroc = auroc / (fpr[0] - fpr[-1])
    auprc = auprc / (tpr[0] - tpr[-1])

    df = pd.DataFrame(
        data = dict(
            cutoff = cutoff,
            tpr = tpr,
            fpr = fpr, 
            precision = precision
        )
    )
    
    res = dict(
        # df : pandas.DataFrame
        #   A DataFrame containing AUROC and AUPRC at various cutoffs.
        df = df,
        
        # AUROC : float
        #   AUROC.
        AUROC = auroc,
        
        # AUPRC : float
        #   AUPRC.
        AUPRC = auprc
    )
    return(res)
