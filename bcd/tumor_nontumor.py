# tumor_nontumor.py - tumor vs. non-tumor classification.

# Inputs
# * CalicoST - spot-wise tumor prop.
# * CopyKAT - label of cell ploidy ('diploid' vs. 'aneuploid');
# * InferCNV - cell x gene CNA expression matrix;
# * Numbat - cell-wise 'tumor' vs. 'normal' assignment;
# * XClone - ...;


import anndata as ad
import functools
import gc
import itertools
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy as sp
import seaborn as sns
import subprocess
import sys
import time
from logging import info, error
from logging import warning as warn
from sccnasim.xlib.xlog import init_logging
from sccnasim.xlib.xrange import format_chrom
from sklearn.cluster import KMeans
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score,
    f1_score, adjusted_rand_score,
    confusion_matrix
)

from .app import APP, VERSION
#APP = "bcd"
#VERSION = "0.4.0"



############################################
#------------------ main ------------------#
############################################

def tumor_nontumor_main(
    sid,
    tool_list,
    out_dir,
    truth_fn,
    tumor_labels,
    overlap_how = 'isec',
    fig_dpi = 300,
    verbose = True
):
    """Main function.

    Parameters
    ----------
    sid : str
        Sample ID.
    tool_list : list of Tool
        A list of tool-specific :class:`~.tool.Tool` objects.
    out_dir : str
        The output folder.
    truth_fn : str
        A header-free file stroing the ground truth.
        Its first two columns should be:
        - `barcode` and `annotation`.
    tumor_labels : str or list of str
        The cell type labels for tumor cells in `truth_fn`.
    overlap_how : {"isec"}
        How to subset the tool matrices given the overlap cells.
        - "isec"
            Subset tool matrix by intersected cells only.
    fig_dpi : int, default 300
        Resolution of the plot.
    verbose : bool, default True
        Whether to show detailed logging information.

    Returns
    -------
    int
        The return code. 0 if success, negative otherwise.
    dict
        Results.
    """
    conf = Config()

    conf.sid = sid
    conf.tool_list = tool_list
    conf.out_dir = out_dir
    conf.truth_fn = truth_fn
    conf.tumor_labels = tumor_labels
    conf.overlap_how = overlap_how

    conf.fig_dpi = fig_dpi
    conf.verbose = verbose

    ret, res = bcd_run(conf)
    return((ret, res))



def bcd_run(conf):
    init_logging(stream = sys.stdout)

    ret = -1
    res = None

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)
    info("%s (VERSION %s)." % (APP, VERSION))

    try:
        res = bcd_core(conf)
    except ValueError as e:
        error(str(e))
        error("Running program failed.")
        error("Quiting ...")
        ret = -1
    else:
        info("All Done!")
        ret = 0
    finally:
        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        info("end time: %s" % time_str)
        info("time spent: %.2fs" % (end_time - start_time, ))

    return((ret, res))



def bcd_core(conf):
    bcd_init(conf)
    info("Configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")


    pp_dir = os.path.join(conf.out_dir, "0_pp")
    os.makedirs(pp_dir, exist_ok = True)


    # predict tumor vs. non-tumor labels.
    info("predict tumor vs. non-tumor labels ...")

    #res_dir = os.path.join(pp_dir, "tools")
    res_dir = pp_dir
    os.makedirs(res_dir, exist_ok = True)
    predict_res = run_predict(
        tool_list = conf.tool_list,
        out_dir = res_dir,
        verbose = conf.verbose
    )


    # extract ground truth.
    info("extract ground truth ...")

    res_dir = os.path.join(pp_dir, "truth")
    os.makedirs(res_dir, exist_ok = True)
    truth_res = run_truth(
        truth_fn = conf.truth_fn,
        out_fn = os.path.join(res_dir, "truth.tsv"),
        tumor_labels = conf.tumor_labels,
        verbose = conf.verbose
    )


    # run core pipeline.
    info("run core pipeline ...")

    res = bcd_core_pipeline(
        sid = conf.sid,
        tool_list = conf.tool_list,
        tool_fn_list = predict_res["out_fns"],
        out_dir = conf.out_dir,
        truth_fn = truth_res['out_fn'],
        overlap_how = conf.overlap_how,
        fig_dpi = conf.fig_dpi,
        verbose = conf.verbose
    )
    return(res)



def bcd_init(conf):
    # check args.
    assert len(conf.tool_list) > 0

    os.makedirs(conf.out_dir, exist_ok = True)
    assert_e(conf.truth_fn)



##############################################
#------------------ config ------------------#
##############################################

class Config:
    def __init__(self):
        self.sid = None
        self.tool_list = None
        self.out_dir = None
        self.truth_fn = None
        self.tumor_labels = None
        self.overlap_how = "isec"
        self.fig_dpi = 300
        self.verbose = True


    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout

        s =  "%s\n" % prefix
        s += "%ssid = %s\n" % (prefix, self.sid)
        s += "%slen(tool_list) = %d\n" % (prefix, len(self.tool_list))
        s += "%stid list = '%s'\n" % (prefix, ", ".join([tool.tid for tool in self.tool_list]))
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%struth_fn = %s\n" % (prefix, self.truth_fn)
        s += "%stumor_labels = %s\n" % (prefix, str(self.tumor_labels))
        s += "%soverlap_how = %s\n" % (prefix, self.overlap_how)
        s += "%sfig_dpi = %d\n" % (prefix, self.fig_dpi)
        s += "%sverbose = %s\n" % (prefix, str(self.verbose))
        s += "%s\n" % prefix

        fp.write(s)



############################################
#------------------ core ------------------#
############################################

def bcd_core_pipeline(
    sid,
    tool_list,
    tool_fn_list,
    out_dir,
    truth_fn,
    overlap_how = 'isec',
    fig_dpi = 300,
    verbose = True
):
    """Pipeline.

    Parameters
    ----------
    sid : str
        Sample ID.
    tool_list : list of Tool
        A list of tool-specific :class:`~.tool.Tool` objects.
    tool_fn_list : list of str
        A list of tool-specific files storing predicted labels.
    out_dir : str
        The output folder.
    truth_fn : str
        An file storing ground truth labels.
    overlap_how
    fig_dpi
    verbose
        See :func:`~.main.bcd_main()` for details.

    Returns
    -------
    dict
        Results.
    """
    # check args.
    if len(tool_list) <= 0:
        info("no input tool data, skip all next steps ...")
        return(dict())

    assert len(tool_list) == len(tool_fn_list)
    for fn in tool_fn_list:
        assert_e(fn)
    os.makedirs(out_dir, exist_ok = True)

    step = 1


    # subset tool and truth labels given overlapping cells.
    info("subset tool and truth labels given overlapping cells ...")

    res_dir = os.path.join(out_dir, "%d_overlap" % step)
    os.makedirs(res_dir, exist_ok = True)
    overlap_res = run_overlap(
        tool_list = tool_list,
        tool_fn_list = tool_fn_list,
        truth_fn = truth_fn,
        overlap_how = overlap_how,
        out_dir = res_dir,
        out_prefix = "overlap",
        verbose = verbose
    )
    step += 1


    # calculate metrics for each tool.
    info("calculate metrics for each tool ...")

    res_dir = os.path.join(out_dir, "%d_metric" % step)
    os.makedirs(res_dir, exist_ok = True)
    metric_res = run_metric(
        tool_list = tool_list,
        tool_fn_list = overlap_res["out_tool_fn_list"],
        truth_fn = overlap_res["out_truth_fn"],
        out_fn = os.path.join(res_dir, "metrics.tsv"),
        verbose = verbose
    )
    step += 1


    # plot metrics.
    info("plot metrics ...")

    res_dir = os.path.join(out_dir, "%d_plot" % step)
    os.makedirs(res_dir, exist_ok = True)
    plot_res = run_plot(
        sid,
        tool_list = tool_list,
        tool_fn_list = overlap_res["out_tool_fn_list"],
        truth_fn = overlap_res["out_truth_fn"],
        metric_fn = metric_res['out_fn'],
        out_dir = res_dir,
        out_prefix = sid,
        fig_dpi = fig_dpi,
        verbose = verbose
    )
    step += 1


    res = plot_res
    return(res)



####################################################
#------------------ steps.metric ------------------#
####################################################

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



# Note that when true labels contain only one class, ARI would be zero
# when test labels has two classes.
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



#####################################################
#------------------ steps.overlap ------------------#
#####################################################

def run_overlap(
    tool_list,
    tool_fn_list,
    truth_fn,
    out_dir,
    out_prefix,
    overlap_how = 'isec',
    verbose = True
):
    """Subset the label files given overlapping cells.

    Parameters
    ----------
    tool_list : list of Tool
        A list of tool-specific :class:`~.tool.Tool` objects.
    tool_fn_list : list of str
        A list of tool-specific label files containing tumor predictions.
    truth_fn : str
        A TSV file storing ground truth of tumor labels.
    out_dir : str
        The output folder.
    out_prefix : str
        The prefix to output files.
    verbose : bool, default True
        Whether to show detailed logging information.

    Returns
    -------
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    return overlap_isec_cells(
        tool_list = tool_list,
        tool_fn_list = tool_fn_list,
        truth_fn = truth_fn,
        out_dir = out_dir,
        out_prefix = out_prefix,
        verbose = verbose
    )



def overlap_isec_cells(
    tool_list,
    tool_fn_list,
    truth_fn,
    out_dir,
    out_prefix,
    verbose = True
):
    """Subset the data by intersected cells only.

    Parameters
    ----------
    tool_list
    tool_fn_list
    truth_fn
    out_dir
    out_prefix
    verbose
        See :func:`run_overlap` for details.

    Returns
    -------
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    # check args.
    assert len(tool_list) > 0
    assert len(tool_fn_list) == len(tool_list)
    for fn in tool_fn_list:
        assert_e(fn)
    assert_e(truth_fn)
    os.makedirs(out_dir, exist_ok = True)


    # get overlap (intersected) cells among tools.
    if verbose:
        info("get overlap (intersected) cells from %d tool files ..." % \
            len(tool_fn_list))

    ovp_cells = None
    for i, (tool, fn) in enumerate(zip(tool_list, tool_fn_list)):
        df = pd.read_csv(fn, delimiter = '\t')
        if i == 0:
            ovp_cells = df["barcode"]
        else:
            ovp_cells = np.intersect1d(df["barcode"], ovp_cells)
    fn = os.path.join(out_dir, "%s.tools.intersect.cells.tsv" % out_prefix)
    np.savetxt(fn, ovp_cells, fmt = "%s", delimiter = '\n')

    if verbose:
        info("tools: %d overlap cells." % len(ovp_cells))


    # further overlap with truth cells.
    if verbose:
        info("further overlap with truth cells ...")

    truth = pd.read_csv(truth_fn, delimiter = '\t')
    ovp_cells = np.intersect1d(truth["barcode"], ovp_cells)
    fn = os.path.join(out_dir,
            "%s.tools_and_truth.intersect.cells.tsv" % out_prefix)
    np.savetxt(fn, ovp_cells, fmt = "%s", delimiter = '\n')

    if verbose:
        info("truth: %d overlap cells." % len(ovp_cells))


    # subset data given overlap cells.
    if verbose:
        info("subset data given overlap cells ...")

    out_tool_fn_list = []
    for tool, tool_fn in zip(tool_list, tool_fn_list):
        tid = tool.tid.lower()
        df = pd.read_csv(tool_fn, delimiter = '\t')
        df.index = df['barcode']
        df = df.loc[ovp_cells].copy()

        fn = os.path.join(out_dir, "%s.%s.tsv" % (out_prefix, tid))
        df.to_csv(fn, sep = '\t', index = False)
        out_tool_fn_list.append(fn)

    truth = pd.read_csv(truth_fn, delimiter = '\t')
    truth.index = truth['barcode']
    truth = truth.loc[ovp_cells].copy()

    out_truth_fn = os.path.join(out_dir, "%s.truth.tsv" % (out_prefix, ))
    truth.to_csv(out_truth_fn, sep = '\t', index = False)


    res = dict(
        # out_tool_fn_list : list of str
        #   A list of output tool-specific label files, in the same order as
        #   `tool_list` and `tool_fn_list`.
        out_tool_fn_list = out_tool_fn_list,

        # out_truth_fn : str
        #   Output subset truth label file.
        out_truth_fn = out_truth_fn
    )
    return(res)



##################################################
#------------------ steps.plot ------------------#
##################################################

def run_plot(
    sid,
    tool_list,
    tool_fn_list,
    truth_fn,
    metric_fn,
    out_dir,
    out_prefix,
    fig_dpi = 300,
    verbose = True
):
    # check args.
    assert len(tool_list) == len(tool_fn_list)
    for fn in tool_fn_list:
        assert_e(fn)
    assert_e(truth_fn)
    os.makedirs(out_dir, exist_ok = True)


    outfig_labels_confusion_matrix = os.path.join(
        out_dir,
        "%s.labels.confusion_matrix.jpg" % (out_prefix,)
    )
    plot_labels_confusion_matrix(
        tool_list = tool_list,
        tool_fn_list = tool_fn_list,
        truth_fn = truth_fn,
        out_fig_fn = outfig_labels_confusion_matrix,
        fig_dpi = fig_dpi,
        verbose = verbose
    )


    outfig_labels_hist = os.path.join(
        out_dir,
        "%s.labels.histogram.jpg" % (out_prefix,)
    )
    plot_labels_hist(
        tool_list = tool_list,
        tool_fn_list = tool_fn_list,
        out_fig_fn = outfig_labels_hist,
        fig_dpi = fig_dpi,
        verbose = verbose
    )


    outfig_metrics_bar = os.path.join(
        out_dir,
        "%s.metrics.bar.jpg" % (out_prefix,)
    )
    plot_metrics_bar(
        metric_fn = metric_fn,
        out_fig_fn = outfig_metrics_bar,
        tid_col = 'tool',
        tid_list = None,
        metric_list = None,
        fig_dpi = fig_dpi,
        verbose = verbose
    )


    outfig_metrics_radar = os.path.join(
        out_dir,
        "%s.metrics.radar.jpg" % (out_prefix,)
    )
    plot_metrics_radar(
        metric_fn = metric_fn,
        out_fig_fn = outfig_metrics_radar,
        fig_dpi = fig_dpi,
        verbose = verbose
    )


    res = dict(
        outfig_labels_confusion_matrix = outfig_labels_confusion_matrix,
        outfig_labels_hist = outfig_labels_hist,
        outfig_metrics_bar = outfig_metrics_bar,
        outfig_metrics_radar = outfig_metrics_radar
    )
    return(res)



def plot_labels_confusion_matrix(
    tool_list,
    tool_fn_list,
    truth_fn,
    out_fig_fn,
    fig_width = None,
    fig_height = None,
    fig_nrow = None,
    fig_ncol = 3,
    fig_dpi = 300,
    verbose = True
):
    # check args.
    assert len(tool_list) == len(tool_fn_list)
    for fn in tool_fn_list:
        assert_e(fn)
    assert_e(truth_fn)


    # plot confusion matrix for each tool.
    if verbose:
        info("plot confusion matrix for each tool's labels ...")

    n = len(tool_list)
    fig_ncol = min(n, fig_ncol)
    if fig_nrow is None:
        fig_nrow = int(np.ceil(n / fig_ncol))
    if fig_width is None:
        fig_width = 2.5 * fig_ncol
    if fig_height is None:
        fig_height = 2.5 * fig_nrow
    fig, axes = plt.subplots(
        fig_nrow, fig_ncol,
        figsize = (fig_width, fig_height)
    )
    if n == 1:
        axes = [axes]

    df = pd.read_csv(truth_fn, sep = '\t')
    truth_labels = df['annotation'].to_numpy()
    for ax, tool, tool_fn in zip(axes, tool_list, tool_fn_list):
        tid = tool.tid
        if verbose:
            info("process %s ..." % tid)

        df = pd.read_csv(tool_fn, sep = '\t')
        tool_labels = df['prediction'].to_numpy()
        cm = confusion_matrix(
            truth_labels, tool_labels,
            labels = ['normal', 'tumor']
        )
        sns.heatmap(
            cm, annot = True, fmt = 'd', cmap = 'Blues',
            xticklabels = ['normal', 'tumor'],
            yticklabels = ['normal', 'tumor'],
            cbar = False,
            ax = ax
        )
        ax.set_title(f'{tid}')
        ax.set_xlabel('Predicted')
        ax.set_ylabel('True')

    plt.tight_layout()
    fig.savefig(out_fig_fn, dpi = fig_dpi, bbox_inches = 'tight')
    plt.close(fig)


    res = dict(
        out_fig_fn = out_fig_fn
    )
    return(res)



def plot_labels_hist(
    tool_list,
    tool_fn_list,
    out_fig_fn,
    fig_width = 3,
    fig_height = 3,
    fig_dpi = 300,
    verbose = True
):
    # check args.
    assert len(tool_list) == len(tool_fn_list)
    for fn in tool_fn_list:
        assert_e(fn)


    # stat.
    tumor_counts, normal_counts = [], []
    for tool, tool_fn in zip(tool_list, tool_fn_list):
        df = pd.read_csv(tool_fn, sep = '\t')
        t = (df['prediction'] == 'tumor').sum()
        tumor_counts.append(t)
        normal_counts.append(df.shape[0] - t)


    # plot histogram for each tool.
    if verbose:
        info("plot histogram for each tool's labels ...")

    fig, ax = plt.subplots(figsize = (fig_width, fig_height))
    bar_width = 0.35
    pos = np.arange(len(tool_list))

    ax.bar(pos - bar_width/2, normal_counts, bar_width,
            label = 'normal', color = '#1f77b4')
    ax.bar(pos + bar_width/2, tumor_counts, bar_width,
            label = 'tumor', color = '#ff7f0e')

    ax.set_xticks(pos)
    ax.set_xticklabels([tool.tid for tool in tool_list])
    ax.set_ylabel('Count')
    ax.set_title('Predicted Label Distribution')
    ax.legend()


    # percentages on top of bars
    total = len(df)
    fongsize = 9
    for i, (n, t) in enumerate(zip(normal_counts, tumor_counts)):
        ax.text(i-bar_width/2, n+total*0.01, f'{n/total:.1%}',
                ha = 'center', va = 'bottom', fontsize = fongsize)
        ax.text(i+bar_width/2, t+total*0.01, f'{t/total:.1%}',
                ha = 'center', va = 'bottom', fontsize = fongsize)

    fig.savefig(out_fig_fn, dpi = fig_dpi, bbox_inches = 'tight')
    plt.close(fig)

    res = dict(
        out_fig_fn = out_fig_fn
    )
    return(res)



def plot_metrics_bar(
    metric_fn,
    out_fig_fn,
    tid_col = 'tool',
    tid_list = None,
    metric_col = 'metric',
    metric_list = None,
    value_col = 'value',
    title = None,
    fig_width = 8,
    fig_height = 8,
    fig_nrow = None,
    fig_ncol = 3,
    fig_dpi = 300,
    verbose = True
):
    # check args.
    df = pd.read_csv(metric_fn, sep = '\t')
    assert tid_col in df.columns
    all_tools = np_unique_keep_order(df[tid_col].to_numpy())
    assert metric_col in df.columns
    all_metrics = np_unique_keep_order(df[metric_col].to_numpy())
    assert value_col in df.columns

    if tid_list is None:
        tid_list = all_tools
    else:
        for tid in tid_list:
            assert tid in all_tools
    if metric_list is None:
        metric_list = all_metrics
    else:
        for metric in metric_list:
            assert metric in all_metrics


    # plot bar for each metric.
    if verbose:
        info("plot bar for each metric ...")

    n = len(metric_list)
    fig_ncol = min(n, fig_ncol)
    if fig_nrow is None:
        fig_nrow = int(np.ceil(n / fig_ncol))
    if fig_width is None:
        fig_width = 2.5 * fig_ncol
    if fig_height is None:
        fig_height = 2.5 * fig_nrow
    fig, axes = plt.subplots(
        fig_nrow, fig_ncol,
        figsize = (fig_width, fig_height),
        sharey = True
    )
    if n == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    bar_width = 0.6
    pos       = np.arange(len(tid_list))

    for ax, metric in zip(axes, metric_list):
        values = []
        for tid in tid_list:
            v = df.loc[(df[tid_col] == tid) & (df[metric_col] == metric), value_col]
            v = v.to_numpy()
            assert len(v) == 1
            values.append(v[0])
        bars = ax.bar(pos, values, bar_width,
                      color = '#1f77b4', edgecolor = 'black')
        ax.set_title(metric, fontsize = 12)
        ax.set_xticks(pos)
        ax.set_xticklabels(tid_list, rotation = 45, ha = 'right')
        ax.set_ylim(0, 1)

        # show value on top of each bar
        for bar, val in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.02,
                f'{val:.3f}',
                ha = 'center', va = 'bottom', fontsize = 9
            )

    fig.supylabel('Score', fontsize = 12)
    plt.tight_layout(rect = [0, 0, 1, 0.95])
    fig.savefig(out_fig_fn, dpi = fig_dpi, bbox_inches = 'tight')
    plt.close(fig)

    res = dict(
        out_fig_fn = out_fig_fn
    )
    return(res)



def plot_metrics_radar(
    metric_fn,
    out_fig_fn,
    tid_col = 'tool',
    tid_list = None,
    metric_col = 'metric',
    metric_list = None,
    value_col = 'value',
    title = None,
    fig_width = 8,
    fig_height = 8,
    fig_dpi = 300,
    linewidth = 2,
    alpha = 0.1,
    verbose = True
):
    # check args.
    df = pd.read_csv(metric_fn, sep = '\t')
    assert tid_col in df.columns
    all_tools = np_unique_keep_order(df[tid_col].to_numpy())
    assert metric_col in df.columns
    all_metrics = np_unique_keep_order(df[metric_col].to_numpy())
    assert value_col in df.columns

    if tid_list is None:
        tid_list = all_tools
    else:
        for tid in tid_list:
            assert tid in all_tools
    if metric_list is None:
        metric_list = all_metrics
    else:
        for metric in metric_list:
            assert metric in all_metrics


    # plot radar for each tool.
    if verbose:
        info("plot radar for each tool's metrics ...")

    angles = np.linspace(
        0, 2*np.pi, len(metric_list),
        endpoint = False
    ).tolist()
    angles += angles[:1]      # close circle

    fig, ax = plt.subplots(
        figsize = (fig_width, fig_height),
        subplot_kw = dict(polar = True)
    )
    for tid in tid_list:
        if verbose:
            info("process %s ..." % tid)

        values = []
        for m in metric_list:
            v = df.loc[(df[tid_col] == tid) & (df[metric_col] == m), value_col]
            v = v.to_numpy()
            assert len(v) == 1
            values.append(v[0])
        values += values[:1]
        ax.plot(angles, values, 'o-', linewidth = linewidth, label = tid)
        ax.fill(angles, values, alpha = alpha)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metric_list)
    ax.set_ylim(0, 1)
    if title is not None:
        ax.set_title(title, size = 16, pad = 20)
    ax.legend(loc = 'upper right', bbox_to_anchor = (1.3, 1.0))

    plt.tight_layout()
    fig.savefig(out_fig_fn, dpi = fig_dpi, bbox_inches = 'tight')
    plt.close(fig)


    res = dict(
        out_fig_fn = out_fig_fn
    )
    return(res)



#####################################################
#------------------ steps.predict ------------------#
#####################################################

def run_predict(
    tool_list, out_dir,
    verbose = True
):
    # check args.
    assert len(tool_list) > 0
    os.makedirs(out_dir, exist_ok = True)

    out_fn_list = []
    for tool in tool_list:
        tid = tool.tid.lower()
        info("predict tumor vs. normal for '%s' ..." % tid)

        res_dir = os.path.join(out_dir, tid)
        os.makedirs(res_dir, exist_ok = True)
        out_fn = os.path.join(res_dir, "%s_predictions.tsv" % tid)

        if tid == "calicost":
            tool.predict(
                out_fn = out_fn,
                verbose = verbose
            )

        elif tid == "copykat":
            tool.predict(
                out_fn = out_fn,
                verbose = verbose
            )

        elif tid == "infercnv":
            out_fn = tool.predict(
                res_dir,
                ref_expr = 1,
                dist = 'euclidean',
                hclust = 'ward.D2',
                cna_score_how = 'mad',
                verbose = verbose
            )

        elif tid == "numbat":
            tool.predict(
                out_fn,
                verbose = verbose
            )

        elif tid == "xclone":
            tool.predict(
                out_fn = out_fn,
                verbose = verbose
            )

        elif tid == "xclone_rdr":
            tool.predict(
                out_fn,
                random_state = 123,
                verbose = verbose
            )
            
        else:
            raise ValueError(f"Error: unknown tool id '{tid}'.")

        out_fn_list.append(out_fn)


    res = dict(
        # out_fns : list of str
        #   Output tumor prediction files.
        #   Each value is a tumor prediction TSV file, in the same order
        #   with `tool_list`.
        out_fns = out_fn_list
    )

    return(res)



###################################################
#------------------ steps.truth ------------------#
###################################################

def run_truth(
    truth_fn,
    out_fn,
    tumor_labels,
    verbose = True
):
    """Format ground truth.

    Read a TSV file containing ground truth annotations, extract barcode
    and true label columns, and save to a TSV file.

    Saves a TSV file with columns: `barcode`, `annotation`.

    Parameters:
    -----------
    truth_fn : str
        The 2-column header-free TSV file with ground truth annotations:
        `barcode` and `annotation`.
    out_fn : str
        File to save the formatted ground truth.
    tumor_labels : str or list of str
        The cell type labels for tumor cells.
    """
    barcode_col = 'barcode'
    anno_col = 'annotation'

    # Check args.
    if not os.path.exists(truth_fn):
        raise ValueError(f"TSV file not found at {truth_fn}!")

    if isinstance(tumor_labels, str):
        tumor_labels = [tumor_labels]

    df = pd.read_csv(truth_fn, delimiter = '\t', header = None)
    df.columns = [barcode_col, anno_col]


    # Create output DataFrame
    out_df = pd.DataFrame({
        barcode_col: df[barcode_col]
    })
    out_df[anno_col] = 'normal'
    out_df.loc[df[anno_col].isin(tumor_labels), anno_col] = 'tumor'


    # Save to TSV
    out_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Ground truth saved to {out_fn}.")


    # Print summary
    n_cells = len(out_df)
    n_unique_labels = len(out_df[anno_col].unique())
    info(f"Extracted {n_cells} cells with ground truth annotations.")
    info(f"Unique labels in {anno_col}: {out_df[anno_col].unique()}.")


    res = dict(
        # out_fn : str
        #   Tumor label ground truth, TSV file.
        out_fn = out_fn
    )
    return(res)



##################################################
#------------------ tools.base ------------------#
##################################################

class Tool:
    def __init__(self, tid):
        self.tid = tid



######################################################
#------------------ tools.calicost ------------------#
######################################################

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
        return calicost_predict_tumor_from_prop(
            tumor_prop_fn = self.tumor_prop_fn,
            out_fn = out_fn,
            prop_col = 'tumor_proportion',
            delimiter = '\t',
            random_state = 123,
            verbose = verbose
        )



def calicost_predict_tumor_from_prop(
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



#####################################################
#------------------ tools.copykat ------------------#
#####################################################

class CopyKAT(Tool):
    def __init__(self, ploidy_pred_fn):
        """CopyKAT object.

        Parameters
        ----------
        ploidy_pred_fn : str
            Path to TSV file with columns 'cell.names' and 'copykat.pred'.
        """
        super().__init__(tid = "CopyKAT")
        self.ploidy_pred_fn = ploidy_pred_fn


    def predict(self, out_fn, verbose = False):
        """Process CopyKAT predictions of tumor vs. non-tumor.

        Read CopyKAT TSV file, rename columns to 'barcode' and 'prediction',
        convert 'diploid' to 'normal' and 'aneuploid' to 'tumor',
        and save to a new TSV file.

        Saves a TSV file with columns:
        - ``barcode``, ``prediction`` ('normal' or 'tumor').
        """
        return copykat_extract_tumor_prediction(
            ploidy_pred_fn = self.ploidy_pred_fn,
            out_fn = out_fn,
            delimiter = '\t',
            verbose = verbose
        )



def copykat_extract_tumor_prediction(
    ploidy_pred_fn,
    out_fn,
    delimiter = '\t',
    verbose = False
):
    # Check args and load data.
    assert_e(ploidy_pred_fn)

    df = pd.read_csv(ploidy_pred_fn, delimiter = delimiter)
    required_cols = ['cell.names', 'copykat.pred']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")

    # sometimes value 'not.defined' exists in column `copykat.pred`.
    valid_preds = {'diploid', 'aneuploid'}
    if not set(df['copykat.pred']).issubset(valid_preds):
        invalid_preds = set(df['copykat.pred']) - valid_preds
        warn(f"Invalid prediction values found: {invalid_preds}."  \
             f"Expected: {valid_preds}")
        df = df.loc[df['copykat.pred'].isin(valid_preds), :].copy()
        warn("%d cells left after removing invalid predictions!" % \
             df.shape[0])


    # Create output DataFrame
    result_df = pd.DataFrame({
        'barcode': df['cell.names'],
        'prediction': df['copykat.pred'].replace(
            {'diploid': 'normal', 'aneuploid': 'tumor'}),
        'ploidy_pred': df['copykat.pred']
    })


    # Save to TSV
    result_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Processed predictions saved to '{out_fn}'.")


    # Print summary
    n_cells = len(df)
    n_tumor = sum(result_df['prediction'] == 'tumor')
    info(f"Processed {n_cells} cells.")
    info(f"Number of tumor cells: {n_tumor}")
    info(f"Number of normal cells: {n_cells - n_tumor}")

    return(out_fn)



######################################################
#------------------ tools.infercnv ------------------#
######################################################

class InferCNV(Tool):
    def __init__(self, obj_fn):
        """InferCNV object.

        Parameters
        ----------
        obj_fn : str
            File storing the inferCNV object.
            Typically using the "MCMC_inferCNV_obj.rds".
        """
        super().__init__(tid = "inferCNV")
        self.obj_fn = obj_fn


    def predict(
        self,
        out_dir,
        ref_expr = 1,
        dist = 'euclidean',
        hclust = 'ward.D2',
        cna_score_how = 'mad',
        verbose = False
    ):
        os.makedirs(out_dir, exist_ok = True)
        out_fn = os.path.join(out_dir, "%s_predictions.tsv" % self.tid.lower())
        res = infercnv_predict_tumor_from_expression(
            obj_fn = self.obj_fn,
            out_fn = out_fn,
            tmp_dir = os.path.join(out_dir, "r2py"),
            ref_expr = ref_expr,
            dist = 'euclidean',
            hclust = 'ward.D2',
            cna_score_how = cna_score_how,
            verbose = verbose
        )
        return out_fn



def infercnv_predict_tumor_from_expression(
    obj_fn,
    out_fn,
    tmp_dir,
    ref_expr = 1,
    dist = 'euclidean',
    hclust = 'ward.D2',
    cna_score_how = 'mad',
    verbose = False
):
    """Predict tumor cells from inferCNV output.

    Read AnnData, perform tumor classification using hierarchical
    clustering on adata.X, and save into TSV file:
    - columns `barcode`, `prediction`.

    Saves TSV files: infercnv_predictions.tsv.
    """
    # Check args.
    assert_e(obj_fn)
    os.makedirs(tmp_dir, exist_ok = True)


    # convert rds to adata
    adata_fn = os.path.join(tmp_dir, 'r2py.h5ad')
    infercnv_extract_cna_expression(obj_fn, adata_fn, tmp_dir = tmp_dir)


    # Perform tumor classification
    # Cluster cells based on expression matrix using hierarchical clustering
    # with (by default) Ward linkage and Euclidean distance,
    # forcing 2 clusters.
    # Label the cluster with the lowest aberration score
    # (mean absolute deviation from median expression) as 'normal' and
    # the other as 'tumor'.
    k = 2
    s  = ""
    s += '''# Classify tumor cells from expression.\n'''
    s += '''\n'''
    s += '''k <- %d\n''' % k
    s += '''ref_expr <- %d\n''' % ref_expr
    s += '''\n'''
    s += '''obj <- readRDS("%s")\n''' % obj_fn
    s += '''mtx <- obj@expr.data\n'''
    s += '''mtx <- t(mtx)         # cell x gene matrix\n'''
    s += '''hc <- hclust(dist(mtx, method = "%s"), method = "%s")\n''' % (dist, hclust)
    s += '''label <- cutree(tree = hc, k = k)\n'''
    s += '''\n'''

    if cna_score_how == 'mad':
        s += '''scores <- rowMeans(abs(mtx - ref_expr))\n'''
    elif cna_score_how == 'md':
        s += '''scores <- rowMeans(mtx - ref_expr)\n'''
    else:
        raise ValueError("unknown cna_score_how '%s'!" % cna_score_how)

    s += '''mean_scores <- rep(Inf, k)\n'''
    s += '''for (i in 1:k)\n'''
    s += '''    mean_scores[i] <- mean(scores[label == i])\n'''
    s += '''normal_cluster <- which.min(mean_scores)\n'''
    s += '''tumor_pred <- ifelse(label == normal_cluster, 'normal', 'tumor')\n'''
    s += '''\n'''
    s += '''df <- data.frame(\n'''
    s += '''    barcode = gsub(".", "-", names(label), fixed = TRUE),\n'''
    s += '''    prediction = tumor_pred,\n'''
    s += '''    cna_score = scores\n'''
    s += ''')\n'''
    s += '''write.table(\n'''
    s += '''    df,\n'''
    s += '''    file = "%s",\n''' % out_fn
    s += '''    sep = "\\t",\n'''
    s += '''    row.names = FALSE,\n'''
    s += '''    col.names = TRUE\n'''
    s += ''')\n'''
    s += '''\n'''

    script_fn = os.path.join(tmp_dir, "classify_tumor_from_expression.R")
    with open(script_fn, "w") as fp:
        fp.write(s)


    # run the R script.
    if verbose:
        info("run the R script to classify tumor cells from expression ...")
    exe_cmdline("Rscript %s" % script_fn)


    # Save predictions TSV
    info(f"Predictions saved to {out_fn}.")


    # Print Summary.
    df = pd.read_csv(out_fn, sep = '\t')
    n_cells = df.shape[0]
    n_tumor = np.sum(df['prediction'] == 'tumor')
    n_normal = np.sum(df['prediction'] == 'normal')
    info("Number of all_cells=%d; tumor_cells=%d; normal_cells=%d." % \
        (n_cells, n_tumor, n_normal))

    return(out_fn)



def infercnv_extract_cna_expression(obj_fn, out_fn, tmp_dir, verbose = False):
    """Extract inferCNV expression matrix and convert it to python object.

    Parameters
    ----------
    out_fn : str
        Output ".h5ad" file storing the cell x gene matrix.
    tmp_dir : str
        The folder to store temporary data.
    verbose : bool, default False
        Whether to show detailed logging information.

    Returns
    -------
    str
        Converted adata file.
    """
    # check args.
    if verbose:
        info("check args ...")
    assert_e(obj_fn)
    os.makedirs(tmp_dir, exist_ok = True)


    # generate R scripts.
    if verbose:
        info("generate R scripts ...")
    cell_fn = os.path.join(tmp_dir, "barcodes.tsv")
    gene_fn = os.path.join(tmp_dir, "genes.tsv")
    mtx_fn = os.path.join(tmp_dir, "matrix.mtx")

    s  = ""
    s += '''# extract cell x gene expression matrix from infercnv output.\n'''
    s += '''\n'''
    s += '''obj <- readRDS("%s")\n''' % obj_fn
    s += '''mtx <- obj@expr.data\n'''
    s += '''mtx <- t(mtx)         # cell x gene matrix\n'''
    s += '''\n'''
    s += '''write(\n'''
    s += '''    rownames(mtx),\n'''
    s += '''    file = "%s"\n''' % cell_fn
    s += ''')\n'''
    s += '''\n'''
    s += '''write(\n'''
    s += '''    colnames(mtx),\n'''
    s += '''    file = "%s"\n''' % gene_fn
    s += ''')\n'''
    s += '''\n'''
    s += '''write.table(\n'''
    s += '''    mtx,\n'''
    s += '''    file = "%s",\n''' % mtx_fn
    s += '''    row.names = FALSE,\n'''
    s += '''    col.names = FALSE\n'''
    s += ''')\n'''
    s += '''\n'''

    script_fn = os.path.join(tmp_dir, "extract_infercnv.R")
    with open(script_fn, "w") as fp:
        fp.write(s)


    # run the R script to save expression matrix into file.
    if verbose:
        info("run the R script to save expression matrix into file ...")
    exe_cmdline("Rscript %s" % script_fn)


    # load matrix into anndata.
    if verbose:
        info("load matrix into anndata and save into .h5ad file ...")
    barcodes = pd.read_csv(cell_fn, header = None)
    barcodes.columns = ["cell"]
    genes = pd.read_csv(gene_fn, header = None)
    genes.columns = ["gene"]

    mtx = np.loadtxt(mtx_fn)
    adata = ad.AnnData(
        X = mtx,
        obs = barcodes,
        var = genes
    )
    save_h5ad(adata, out_fn)
    if verbose:
        info("saved adata shape = %s." % str(adata.shape))

    return(out_fn)



####################################################
#------------------ tools.numbat ------------------#
####################################################

class Numbat(Tool):
    def __init__(self, clone_post_fn):
        """Numbat object.

        Parameters
        ----------
        clone_post_fn : str
            Path to Numbat TSV file with columns including 'cell' and
            'compartment_opt'.
        """
        super().__init__(tid = "Numbat")
        self.clone_post_fn = clone_post_fn


    def predict(self, out_fn, verbose = False):
        """Predict tumor cells from Numbat output.

        Saves a TSV file with columns: `barcode`, `prediction`
        ('normal' or 'tumor') to out_dir/numbat_predictions.tsv.
        """
        return numbat_extract_tumor(
            clone_post_fn = self.clone_post_fn,
            out_fn = out_fn,
            barcode_col = 'cell',
            label_col = 'compartment_opt',
            delimiter = '\t',
            verbose = verbose
        )



def numbat_extract_tumor(
    clone_post_fn,
    out_fn,
    barcode_col = 'cell',
    label_col = 'compartment_opt',
    p_cnv_col = 'p_cnv',
    delimiter = '\t',
    verbose = False
):
    # Check args.
    assert_e(clone_post_fn)

    df = pd.read_csv(clone_post_fn, delimiter = delimiter)
    required_cols = [barcode_col, label_col, p_cnv_col]
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")


    # Filter out rows with empty or invalid labels.
    n_cells_init = len(df)
    df = df.dropna(subset = [label_col])
    df = df.loc[df[label_col].isin(['normal', 'tumor'])].copy()
    n_cells = len(df)
    n_removed = n_cells_init - n_cells
    if n_removed > 0:
        warn(f"Removed %d cells with empty/invalid '%s' values." % \
            (n_removed, label_col))
    if n_cells == 0:
        error(f"No cells remain after removing empty/invalid '%s' values." % \
              label_col)
        raise ValueError


    # Save to TSV
    result_df = pd.DataFrame({
        'barcode': df[barcode_col].to_numpy(),
        'prediction': df[label_col].to_numpy(),
        'p_cnv': df[p_cnv_col].to_numpy()
    })
    result_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Predictions saved to {out_fn}.")


    # Print Summary.
    df = pd.read_csv(out_fn, sep = '\t')
    n_cells = df.shape[0]
    n_tumor = np.sum(df['prediction'] == 'tumor')
    n_normal = np.sum(df['prediction'] == 'normal')
    info("Number of all_cells=%d; tumor_cells=%d; normal_cells=%d." % \
        (n_cells, n_tumor, n_normal))

    return(out_fn)



####################################################
#------------------ tools.xclone ------------------#
####################################################

class XClone(Tool):
    def __init__(self, xclone_tumor_pred_fn):
        """XClone object.

        Parameters
        ----------
        xclone_tumor_pred_fn : str
            Path to XClone TSV file with columns 'barcode' and 'prediction'.
        """
        super().__init__(tid = "XClone")
        self.xclone_tumor_pred_fn = xclone_tumor_pred_fn


    def predict(self, out_fn, verbose = False):
        """Extract tumor predictions from XClone output.

        Reads XClone TSV file with barcode and prediction columns,
        validates the data, and saves to a new TSV file.

        Saves a TSV file with columns:
        - ``barcode``, ``prediction`` ('normal' or 'tumor').
        """
        return xclone_extract_tumor_prediction(
            xclone_tumor_pred_fn = self.xclone_tumor_pred_fn,
            out_fn = out_fn,
            delimiter = '\t',
            verbose = verbose
        )



def xclone_extract_tumor_prediction(
    xclone_tumor_pred_fn,
    out_fn,
    delimiter = '\t',
    verbose = False
):
    """Extract tumor predictions from XClone output.

    Parameters
    ----------
    xclone_tumor_pred_fn : str
        Path to XClone TSV file containing columns: 'barcode' and 'prediction'.
    out_fn : str
        Output file path for tumor predictions.
    delimiter : str, default '\\t'
        Delimiter used in the input TSV file.
    verbose : bool, default False
        Whether to show detailed logging information.

    Returns
    -------
    str
        Path to output prediction file.
    """
    # Check args and load data.
    assert_e(xclone_tumor_pred_fn)

    df = pd.read_csv(xclone_tumor_pred_fn, delimiter = delimiter)
    required_cols = ['barcode', 'prediction']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")

    # Validate prediction values
    valid_preds = {'normal', 'tumor'}
    if not set(df['prediction']).issubset(valid_preds):
        invalid_preds = set(df['prediction']) - valid_preds
        warn(f"Invalid prediction values found: {invalid_preds}. "
             f"Expected: {valid_preds}")
        df = df.loc[df['prediction'].isin(valid_preds), :].copy()
        warn(f"{len(df)} cells left after removing invalid predictions!")

    # Create output DataFrame
    result_df = pd.DataFrame({
        'barcode': df['barcode'],
        'prediction': df['prediction']
    })


    # Save to TSV
    result_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Predictions saved to '{out_fn}'.")


    # Print summary
    n_cells = len(result_df)
    n_tumor = sum(result_df['prediction'] == 'tumor')
    n_normal = n_cells - n_tumor
    info(f"Processed {n_cells} cells.")
    info(f"Number of tumor cells: {n_tumor}")
    info(f"Number of normal cells: {n_normal}")

    return(out_fn)



##################################################
#------------------ utils.base ------------------#
##################################################

def assert_e(path):
    """Assert file or folder exists, mimicking shell "test -e"."""
    assert path is not None and os.path.exists(path)



def exe_cmdline(cmdline):
    """Execulate command line."""
    proc = subprocess.Popen(
        args = cmdline,
        shell = True,
        executable = "/bin/bash",
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    outs, errs = proc.communicate()
    ret = proc.returncode
    if ret != 0:
        raise RuntimeError(str(errs.decode()))



def expand_grid(d):
    """An equivalent to R expand.grid().

    Parameters
    ----------
    d : dict
        Key of dict is the factor name (column name in returned DataFrame),
        value of dict is a list/vector of factor-specific values.

    Returns
    -------
    pandas.DataFrame
        The combination of all factors.
    """
    # ref: https://stackoverflow.com/questions/71116437/expand-grid-equivalent-to-get-pandas-data-frame-for-prediction-in-python/71376414
    try:
        res = pd.MultiIndex.from_product(
            d.values(),
            names = d.keys()
        ).to_frame().reset_index(drop = True)
        return(res)
    except:
        res = pd.DataFrame(itertools.product(*d.values()), columns = d.keys())
        return(res)


def np_unique_keep_order(x):
    b, idx = np.unique(x, return_index = True)
    return b[np.argsort(idx)]



def is_file_empty(fn):
    """Test whether file is empty."""
    assert os.path.exists(fn)
    return(os.path.getsize(fn) <= 0)



################################################
#------------------ utils.io ------------------#
################################################

def load_h5ad(fn):
    return(ad.read_h5ad(fn))

def save_h5ad(adata, fn):
    return(adata.write_h5ad(fn, compression = "gzip"))
