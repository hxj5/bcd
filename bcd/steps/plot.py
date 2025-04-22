# plot.py


import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from logging import info
from ..utils.base import assert_e



def run_plot(
    sid,
    cna_type,
    out_dir,
    out_prefix,
    roc_fn,
    auroc_fn,
    prc_fn,
    auprc_fn,
    fig_width = 4.25,
    fig_height = 3.25,
    fig_dpi = 300,
    fig_dec = 3,
    fig_legend_xmin = 0.5,
    fig_legend_ymin = 0.12,
    verbose = True
):
    """Main function of plotting ROC and PRC.
    
    Parameters
    ----------
    sid : str
        Sample ID.
    cna_type : str
        CNA types, one of {"gain", "loss", "loh"}.
    out_dir : str
        The output folder.
    out_prefix : str
        The prefix to output files.
    roc_fn : str
        File storing DataFrame that contains ROC at various cutoffs.
        It has four columns "tool", "cutoff", tpr", "fpr".
    auroc_fn : str
        File storing DataFrame that contains AUROC of each tool.
        It has two columns "tool" and "AUROC".
    prc_fn : str
        File storing DataFrame that contains PRC at various cutoffs.
        It has four columns "tool", "cutoff", "precision", "recall".
    auprc_fn : str
        File storing DataFrame that contains AUPRC of each tool.
        It has two columns "tool" and "AUPRC".
    fig_width
    fig_height
    fig_dpi
    fig_dec
    fig_legend_xmin
    fig_legend_ymin
    verbose
        See :func:`~.main.bcd_main()` for details.
        
    Returns
    -------
    Void.
    """
    # check args.
    assert cna_type in ("gain", "loss", "loh")
    os.makedirs(out_dir, exist_ok = True)
    assert_e(roc_fn)
    assert_e(auroc_fn)
    assert_e(prc_fn)
    assert_e(auprc_fn)
    
    cna_name_map = {
        "gain" : "Copy Gain",
        "loss" : "Copy Loss",
        "loh" : "LOH"
    }
        
    out_fns = dict()
    for metric, df_fn, auc_fn in zip(("roc", "prc"), (roc_fn, prc_fn), \
                                     (auroc_fn, auprc_fn)):
        if verbose:
            info("plot %s ..." % metric.upper())

        out_fn = os.path.join(out_dir, "%s.%s.%s.jpg" % (out_prefix, cna_type, metric))
        plot_metric(
            metric = metric,
            df_fn = df_fn,
            auc_fn = auc_fn,
            out_fn = out_fn,
            title = "%s %s Curve for %s" % (sid, metric.upper(), cna_name_map[cna_type]),
            fig_width = fig_width,
            fig_height = fig_height,
            fig_dpi = fig_dpi,
            fig_dec = fig_dec, 
            fig_legend_xmin = fig_legend_xmin,
            fig_legend_ymin = fig_legend_ymin
        )
        out_fns[metric] = out_fn
        

    res = dict(
        # out_fns : dict of {str : str}
        #   Output figures for each metric.
        #   Kyes are metrics, values are output figure files.
        out_fns = out_fns
    )
    return(res)



def plot_metric(
    metric,
    df_fn,
    auc_fn,
    out_fn,
    title = None,
    fig_width = 4.25,
    fig_height = 3.25,
    fig_dpi = 300,
    fig_dec = 3,
    fig_legend_xmin = 0.5,
    fig_legend_ymin = 0.12,
):
    """Plot ROC or PRC curve.
    
    Parameters
    ----------
    metric : {"roc", "prc"}
        Evaluation metric.
    df_fn : str
        File storing DataFrame that contains ROC or PRC at various cutoffs.
        If ROC, it has four columns "tool", "cutoff", tpr", "fpr";
        otherwise, it has four columns "tool", "cutoff", "precision", "recall".
    auc_fn : str
        File storing DataFrame that contains AUROC or AUPRC of each tool.
        If AUROC, it has two columns "tool" and "AUROC";
        otherwise, it has two columns "tool" and "AUPRC".
    out_fn : str
        Output figure file.
    title : str or None
        Title of the plot figure. None means no title.
    fig_width
    fig_height
    fig_dpi
    fig_dec
    fig_legend_xmin
    fig_legend_ymin
        See :func:`~.main.bcd_main()` for details.
        
    Returns
    -------
    dict
        Results.
    """
    # check args.
    assert metric in ("roc", "prc")
    assert_e(df_fn)
    assert_e(auc_fn)
    
    
    df = pd.read_csv(df_fn, sep = "\t")
    auc = pd.read_csv(auc_fn, sep = "\t")
    
    if metric == "roc":
        for c in ("tool", "cutoff", "tpr", "fpr"):
            assert c in df.columns
    else:
        for c in ("tool", "cutoff", "precision", "recall"):
            assert c in df.columns
            
    if metric == "roc":
        for c in ("tool", "AUROC"):
            assert c in auc.columns
    else:
        for c in ("tool", "AUPRC"):
            assert c in auc.columns
            
    assert len(df["tool"].unique()) == len(auc["tool"].unique())
    assert np.all(np.sort(df["tool"].unique()) == np.sort(auc["tool"].unique()))
            
        
    # prepare plot data.
    df = df.merge(auc, on = "tool", how = "inner")
    
    if fig_dec == 4:
        if metric == "roc":
            df["legend"] = df.apply(lambda x: "%s: AUC=%.4f" % \
                                    (x["tool"], x["AUROC"]), axis = 1)
        else:
            df["legend"] = df.apply(lambda x: "%s: AUC=%.4f" % \
                                    (x["tool"], x["AUPRC"]), axis = 1)
    else:
        if metric == "roc":
            df["legend"] = df.apply(lambda x: "%s: AUC=%.3f" % \
                                    (x["tool"], x["AUROC"]), axis = 1)
        else:
            df["legend"] = df.apply(lambda x: "%s: AUC=%.3f" % \
                                    (x["tool"], x["AUPRC"]), axis = 1)

    
    # plot.
    title_fontsize = None
    axis_fontsize = None
    legend_fontsize = 8
    linewidth = 0.6
    
    plt.figure(figsize = (fig_width, fig_height))
    
    if metric == "roc":
        sns.lineplot(data = df, x = 'fpr', y = 'tpr', hue = 'legend', 
                        linewidth = linewidth, estimator = None, sort = False)
        plt.xlabel("False Positive Rate (1 - Specificity)", fontsize = axis_fontsize)
        plt.ylabel("True Positive Rate (Sensitivity)", fontsize = axis_fontsize)
    else:
        sns.lineplot(data = df, x = 'recall', y = 'precision', hue = 'legend', 
                        linewidth = linewidth, estimator = None, sort = False)
        plt.xlabel("Recall", fontsize = axis_fontsize)
        plt.ylabel("Precision", fontsize = axis_fontsize)
    
    if title is not None:
        plt.title(title)

    if metric == "roc":
        plt.legend(
            loc = 'lower left', 
            bbox_to_anchor=(fig_legend_xmin, fig_legend_ymin), 
            fontsize = legend_fontsize
        )
    else:
        plt.legend(
            loc = 'lower right', 
            bbox_to_anchor=(fig_legend_xmin, fig_legend_ymin), 
            fontsize = legend_fontsize
        )

    plt.grid(False)
    plt.tight_layout()
    
    plt.savefig(out_fn, dpi = fig_dpi)

