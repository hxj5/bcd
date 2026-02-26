# subclonal_structure.py - subclonal structure identification.

# Inputs
# * CalicoST - spot-wise clone label.
# * CopyKAT - hclustering results;
# * InferCNV - cell x gene CNA expression matrix;
# * Numbat - cell-wise clone label;
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
from sklearn.metrics import adjusted_rand_score, confusion_matrix

from .app import APP, VERSION
#APP = "bcd"
#VERSION = "0.4.0"



############################################
#------------------ main ------------------#
############################################

def subclonal_structure_main(
    sid,
    tool_list,
    out_dir,
    truth_fn,
    n_cluster,
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
    n_cluster : int
        Number of clusters for tools that do not output clone labels.
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
    conf.n_cluster = n_cluster
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

    
    # predict subclonal structure.
    info("predict subclonal structure ...")
    
    #res_dir = os.path.join(pp_dir, "tools")
    res_dir = pp_dir
    os.makedirs(res_dir, exist_ok = True)
    predict_res = run_predict(
        tool_list = conf.tool_list, 
        out_dir = res_dir,
        k = conf.n_cluster,
        truth_fn = conf.truth_fn,
        verbose = conf.verbose
    )
    
    
    # extract ground truth.
    info("extract ground truth ...")

    res_dir = os.path.join(pp_dir, "truth")
    os.makedirs(res_dir, exist_ok = True)
    truth_res = run_truth(
        truth_fn = conf.truth_fn, 
        out_fn = os.path.join(res_dir, "truth.tsv"),
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
        self.n_cluster = None
        self.overlap_how = "isec"
        self.fig_dpi = 300
        self.verbose = True


    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
            
        s =  "%s\n" % prefix
        s += "%ssid = %s\n" % (prefix, self.sid)
        s += "%slen(tool_list) = %d\n" % (prefix, len(self.tool_list))
        s += "%stid list = '%s'\n" % (prefix, ", ".join([tool.display_name() for tool in self.tool_list]))
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%struth_fn = %s\n" % (prefix, self.truth_fn)
        s += "%sn_cluster = %s\n" % (prefix, str(self.n_cluster))
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

    ari_list = []

    df = pd.read_csv(truth_fn, sep = '\t')
    truth_labels = df['annotation'].to_numpy()
    for i, (tool, tool_fn) in enumerate(zip(tool_list, tool_fn_list)):
        display_name = tool.display_name()
        if verbose:
            info("process %s ..." % display_name)
        
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
            tool = [tool.display_name() for tool in tool_list],
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



#####################################################
#------------------ steps.overlap ------------------#
#####################################################

def _tool_file_id(tool):
    """Get unique filesystem identifier for tool (for output filenames)."""
    if tool.run_id is not None and str(tool.run_id).strip() != "":
        return "%s_%s" % (tool.tid.lower(), str(tool.run_id))
    return tool.tid.lower()


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
        A list of tool-specific label files containing clone predictions.
    truth_fn : str
        A TSV file storing ground truth of clone labels.
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
        file_id = _tool_file_id(tool)
        df = pd.read_csv(tool_fn, delimiter = '\t')
        df.index = df['barcode']
        df = df.loc[ovp_cells].copy()
        
        fn = os.path.join(out_dir, "%s.%s.tsv" % (out_prefix, file_id))
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


    res = dict(
        outfig_labels_confusion_matrix = outfig_labels_confusion_matrix,
        outfig_metrics_bar = outfig_metrics_bar
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
    # Handle different return types from plt.subplots()
    # plt.subplots() can return: single Axes (n=1), 1D array (1 row/col), or 2D array
    if isinstance(axes, np.ndarray):
        axes = axes.flatten()
    else:
        axes = [axes]
    
    # Hide unused subplots if we have more subplots than tools
    for i in range(n, len(axes)):
        axes[i].set_visible(False)

    df = pd.read_csv(truth_fn, sep = '\t')
    truth_labels = df['annotation'].to_numpy()
    for ax, tool, tool_fn in zip(axes[:n], tool_list, tool_fn_list):
        display_name = tool.display_name()
        if verbose:
            info("process %s ..." % display_name)

        df = pd.read_csv(tool_fn, sep = '\t')
        tool_labels = df['prediction'].to_numpy()
        cm = confusion_matrix(
            truth_labels, tool_labels
        )
        sns.heatmap(
            cm, annot = True, fmt = 'd', cmap = 'Blues',
            cbar = False,
            ax = ax
        )
        ax.set_title(f'{display_name}')
        ax.set_xlabel('Predicted')
        ax.set_ylabel('True')

    plt.tight_layout()
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



#####################################################
#------------------ steps.predict ------------------#
#####################################################


def merge_clones_to_k(
    tool_fn,
    k,
    truth_fn = None,
    method = 'truth',
    verbose = False
):
    """If a tool outputs more than k clones, merge clones to produce exactly k.

    When method='truth', merges the pair of clones with highest similarity
    with respect to the ground truth (similar distribution over truth labels).
    When method='size' or truth_fn is missing, falls back to merging the two
    smallest clusters. Preserves barcode and other columns; only 'prediction'
    is updated.

    Parameters
    ----------
    tool_fn : str
        Path to TSV with columns at least: barcode, prediction.
    k : int
        Target number of clones.
    truth_fn : str or None
        Path to ground truth TSV (barcode, annotation). Required for method='truth'.
        Can be header-free (first two columns taken as barcode, annotation).
    method : {"truth", "size"}
        Merge strategy. "truth": merge the two clones with highest similarity
        to the truth (cosine similarity of count vectors over truth labels).
        "size": repeatedly merge the two smallest clusters.
    verbose : bool
        Whether to log merge steps.

    Returns
    -------
    str
        Path to the same file (overwritten).
    """
    assert_e(tool_fn)
    df = pd.read_csv(tool_fn, sep = '\t')
    assert 'prediction' in df.columns
    pred = df['prediction'].to_numpy()
    unique_labels = np.unique(pred)
    n_clus = len(unique_labels)
    if n_clus <= k:
        return tool_fn

    # Resolve method: use truth only if truth_fn is provided
    use_truth = (method == 'truth' and truth_fn is not None and os.path.exists(truth_fn))
    if method == 'truth' and not use_truth:
        if verbose:
            warn("merge_clones_to_k: truth_fn missing or invalid, falling back to method='size'")
        use_truth = False

    if verbose:
        info("merge_clones_to_k: %d clusters -> %d (method=%s)" % (n_clus, k, 'truth' if use_truth else 'size'))

    # Build label -> size (used for size method and for maintaining rep sizes)
    _, counts = np.unique(pred, return_counts = True)
    label_to_size = dict(zip(unique_labels.tolist(), counts.tolist()))

    # Union-find style
    label_to_rep = {l: l for l in unique_labels.tolist()}
    rep_to_labels = {l: {l} for l in unique_labels.tolist()}
    rep_to_size = {l: label_to_size[l] for l in unique_labels.tolist()}

    # For truth method: overlap with truth and compute per-(rep) count vector over truth labels
    rep_to_truth_vec = None
    if use_truth:
        peek = pd.read_csv(truth_fn, sep = '\t', nrows = 1)
        if peek.columns[0] in ('barcode', 'Barcode', 'BARCODE', 'cell', 'Cell'):
            truth_df = pd.read_csv(truth_fn, sep = '\t').iloc[:, :2]
        else:
            truth_df = pd.read_csv(truth_fn, sep = '\t', header = None).iloc[:, :2]
        truth_df.columns = ['barcode', 'annotation']
        overlap = df[['barcode', 'prediction']].merge(
            truth_df, on = 'barcode', how = 'inner'
        )
        if len(overlap) == 0:
            if verbose:
                warn("merge_clones_to_k: no overlap with truth, falling back to method='size'")
            use_truth = False
        else:
            truth_labels = np.unique(overlap['annotation'].to_numpy())
            truth_to_idx = {t: i for i, t in enumerate(truth_labels)}
            # For each predicted label (we'll use rep as key), count vector over truth
            def count_vec(labels_subset):
                sub = overlap[overlap['prediction'].isin(labels_subset)]
                vec = np.zeros(len(truth_labels))
                for ann, cnt in sub['annotation'].value_counts().items():
                    idx = truth_to_idx.get(ann, -1)
                    if idx >= 0:
                        vec[idx] = cnt
                return vec
            rep_to_truth_vec = {r: count_vec(rep_to_labels[r]) for r in rep_to_size}

    def cosine_similarity(u, v):
        nu = np.asarray(u, dtype = float)
        nv = np.asarray(v, dtype = float)
        norm_u = np.sqrt((nu * nu).sum())
        norm_v = np.sqrt((nv * nv).sum())
        if norm_u == 0 or norm_v == 0:
            return 0.0
        return float(np.dot(nu, nv) / (norm_u * norm_v))

    def pair_to_merge():
        if use_truth and rep_to_truth_vec is not None:
            # Merge the two reps with highest cosine similarity of truth count vectors
            reps = list(rep_to_size.keys())
            best_sim = -2.0
            best_pair = (reps[0], reps[1])
            for i in range(len(reps)):
                for j in range(i + 1, len(reps)):
                    ri, rj = reps[i], reps[j]
                    sim = cosine_similarity(rep_to_truth_vec[ri], rep_to_truth_vec[rj])
                    if sim > best_sim:
                        best_sim = sim
                        best_pair = (ri, rj)
            return best_pair
        else:
            # Two smallest
            sorted_reps = sorted(rep_to_size.keys(), key = lambda r: rep_to_size[r])
            return sorted_reps[0], sorted_reps[1]

    while len(rep_to_size) > k:
        rep_a, rep_b = pair_to_merge()
        # Merge rep_a into rep_b
        for l in rep_to_labels[rep_a]:
            label_to_rep[l] = rep_b
        rep_to_labels[rep_b] |= rep_to_labels[rep_a]
        rep_to_size[rep_b] = rep_to_size[rep_b] + rep_to_size[rep_a]
        if use_truth and rep_to_truth_vec is not None:
            rep_to_truth_vec[rep_b] = np.asarray(rep_to_truth_vec[rep_b]) + np.asarray(rep_to_truth_vec[rep_a])
            del rep_to_truth_vec[rep_a]
        del rep_to_labels[rep_a]
        del rep_to_size[rep_a]

    # Map old prediction to representative, then renumber to 0..k-1
    reps = sorted(rep_to_size.keys())
    rep_to_new = {r: i for i, r in enumerate(reps)}
    new_pred = np.array([rep_to_new[label_to_rep[p]] for p in pred])
    df['prediction'] = new_pred
    df.to_csv(tool_fn, sep = '\t', index = False)
    if verbose:
        info("merge_clones_to_k: wrote %s" % tool_fn)
    return tool_fn


def run_predict(
    tool_list, out_dir, 
    k,
    truth_fn = None,
    verbose = True
):
    # check args.
    assert len(tool_list) > 0
    os.makedirs(out_dir, exist_ok = True)
    
    out_fn_list = []
    for tool in tool_list:
        tid = tool.tid.lower()
        file_id = _tool_file_id(tool)
        info("predict subclonal structure for '%s' ..." % tool.display_name())

        res_dir = os.path.join(out_dir, tid)
        os.makedirs(res_dir, exist_ok = True)
        out_fn = os.path.join(res_dir, "%s_predictions.tsv" % file_id)
        
        if tid == "calicost":
            tool.predict(
                out_fn = out_fn,
                verbose = verbose
            )

        elif tid == "copykat":
            tool.predict(
                out_fn = out_fn,
                k = k,
                tmp_dir = res_dir,
                verbose = verbose
            )
        
        elif tid == "infercnv":
            out_fn = tool.predict(
                out_fn = out_fn,
                k = k,
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
        
        else:
            raise ValueError(f"Error: unknown tool id '{tid}'.")

        # If tool produced more than k clones, merge down to k (by similarity to truth)
        merge_clones_to_k(out_fn, k, truth_fn = truth_fn, method = 'truth', verbose = verbose)
        out_fn_list.append(out_fn)


    res = dict(
        # out_fns : list of str
        #   Output subclonal structure prediction files.
        #   Each value is a prediction TSV file, in the same order
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
    """
    barcode_col = 'barcode'
    anno_col = 'annotation'

    # Check args.
    if not os.path.exists(truth_fn):
        raise ValueError(f"TSV file not found at {truth_fn}!")

    df = pd.read_csv(truth_fn, delimiter = '\t', header = None)
    df.columns = [barcode_col, anno_col]


    # Create output DataFrame
    labels_old = np.unique(df[anno_col])
    label_new = np.arange(len(labels_old))
    label_map = {o:n for o, n in zip(labels_old, label_new)}

    out_df = pd.DataFrame({
        barcode_col: df[barcode_col],
        anno_col: df[anno_col].map(label_map),
        anno_col+"_old": df[anno_col]
    })


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
    def __init__(self, tid, run_id = None):
        self.tid = tid
        self.run_id = run_id

    def display_name(self):
        """Return display name for plots/legends. Includes run_id when present."""
        if self.run_id is not None and str(self.run_id).strip() != "":
            return "%s_%s" % (self.tid, self.run_id)
        return self.tid



######################################################
#------------------ tools.calicost ------------------#
######################################################

class CalicoST(Tool):
    def __init__(self, clone_label_fn, run_id = None):
        """CalicoST Object.
        
        Parameters
        ----------
        clone_label_fn : str
            Path to CalicoST TSV file containing columns: 
            ``BARCODES`` and ``clone_label``.
        run_id : str or None, default None
            Optional run identifier for multiple runs of the same tool.
        """
        super().__init__(tid = "CalicoST", run_id = run_id)
        self.clone_label_fn = clone_label_fn


    def predict(self, out_fn, verbose = False):
        """Extract the subclonal labels from CalicoST output.

        Saves a TSV file with columns: `barcode`, `prediction`,
        and `clone_label`.
        """
        return extract_clonal_labels(
            clone_label_fn = self.clone_label_fn,
            out_fn = out_fn,
            delimiter = '\t',
            verbose = verbose
        )



def extract_clonal_labels(
    clone_label_fn,
    out_fn,
    delimiter = '\t',
    verbose = False
):
    # Check args and load data.
    df = pd.read_csv(clone_label_fn, delimiter = delimiter)

    required_cols = ['BARCODES', 'clone_label']
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")


    # Filter out rows with empty clone labels.
    initial_n_cells = len(df)
    df = df.dropna(subset = ['clone_label'])
    n_cells = len(df)
    n_removed = initial_n_cells - n_cells
    if n_removed > 0:
        warn("Removed %d cells with empty 'clone_label' values." %  \
            (n_removed, ))
    if n_cells == 0:
        error("No cells remain after removing empty 'clone_label' values.")
        raise ValueError


    # Format labels.
    labels_old = np.unique(df['clone_label'])
    label_new = np.arange(len(labels_old))
    label_map = {o:n for o, n in zip(labels_old, label_new)}
    

    # Save to TSV
    result_df = pd.DataFrame({
        'barcode': df['BARCODES'],
        'prediction': df['clone_label'].map(label_map),
        'clone_label': df['clone_label']
    })
    result_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Predictions saved to '{out_fn}'.")

    
    # Print summary
    info(f"Processed {n_cells} cells after filtering.")
    info("In total %d clones." % len(label_new))
    
    return(out_fn)



#####################################################
#------------------ tools.copykat ------------------#
#####################################################

class CopyKAT(Tool):
    def __init__(self, hclust_fn, run_id = None):
        """CopyKAT object.
        
        Parameters
        ----------
        hclust_fn : str
            Path to RDS file containing CopyKAT hclustering results.
        run_id : str or None, default None
            Optional run identifier for multiple runs of the same tool.
        """
        super().__init__(tid = "CopyKAT", run_id = run_id)
        self.hclust_fn = hclust_fn

        
    def predict(self, out_fn, k, tmp_dir, verbose = False):
        """Process CopyKAT predictions of subclonal structures.

        Saves a TSV file with columns: 
        - ``barcode``, ``prediction``.
        """
        return predict_subclones_from_copykat_hclust(
            hclust_fn = self.hclust_fn,
            out_fn = out_fn,
            k = k,
            tmp_dir = tmp_dir,
            verbose = verbose
        )
        
        
        
def predict_subclones_from_copykat_hclust(
    hclust_fn,
    out_fn,
    k,
    tmp_dir,
    verbose = False
):
    # Check args and load data.
    assert_e(hclust_fn)
    
    s = ""
    s += '''# Extract the subclonal structure using cutree on hclust.\n'''
    s += '''\n'''
    s += '''obj <- readRDS("%s")\n''' % hclust_fn
    s += '''label <- cutree(tree = obj, k = %d)\n''' % k
    s += '''df <- data.frame(\n'''
    s += '''    barcode = gsub(".", "-", names(label), fixed = TRUE),\n'''
    s += '''    prediction = label - 1\n'''
    s += ''')\n'''
    s += '''write.table(\n'''
    s += '''    df,\n'''
    s += '''    file = "%s",\n''' % out_fn
    s += '''    sep = "\\t",\n'''
    s += '''    row.names = FALSE,\n'''
    s += '''    col.names = TRUE\n'''
    s += ''')\n'''
    s += '''\n'''
    
    script_fn = os.path.join(tmp_dir, "predict_subclones_from_copykat_hclust.R")
    with open(script_fn, "w") as fp:
        fp.write(s)


    # run the R script.
    if verbose:
        info("run the R script to predict subclones from hclust ...")

    exe_cmdline("Rscript %s" % script_fn)

    
    # Save to TSV
    info(f"Processed predictions saved to '{out_fn}'.")


    # Print summary
    df = pd.read_csv(out_fn, sep = '\t')
    labels = np.unique(df['prediction'])
    n_cells = len(df)
    info(f"Processed {n_cells} cells after filtering.")
    info("In total %d clones." % len(labels))
    
    return(out_fn)



######################################################
#------------------ tools.infercnv ------------------#
######################################################

class InferCNV(Tool):
    def __init__(self, obj_fn, run_id = None):
        """InferCNV object.
        
        Parameters
        ----------
        obj_fn : str
            File storing the inferCNV object. 
            Typically using the "MCMC_inferCNV_obj.rds".
        run_id : str or None, default None
            Optional run identifier for multiple runs of the same tool.
        """
        super().__init__(tid = "inferCNV", run_id = run_id)
        self.obj_fn = obj_fn
        
        
    def predict(
        self,
        out_fn,
        k,
        verbose = False
    ):
        out_dir = os.path.dirname(out_fn)
        os.makedirs(out_dir, exist_ok = True)
        res = predict_subclones_from_expression(
            obj_fn = self.obj_fn,
            out_fn = out_fn,
            k = k,
            tmp_dir = os.path.join(out_dir, "r2py"),
            verbose = verbose
        )
        return out_fn



def predict_subclones_from_expression(
    obj_fn,
    out_fn,
    k,
    tmp_dir,
    dist = 'euclidean',
    hclust = 'ward.D2',
    verbose = False
):
    # Check args.
    assert_e(obj_fn)
    os.makedirs(tmp_dir, exist_ok = True)
    
    
    # convert rds to adata
    adata_fn = os.path.join(tmp_dir, 'r2py.h5ad')
    extract_cna_expression(obj_fn, adata_fn, tmp_dir = tmp_dir)
    
    
    # Perform subclone identification.
    s  = ""
    s += '''# Identify subclones from expression.\n'''
    s += '''\n'''
    s += '''obj <- readRDS("%s")\n''' % obj_fn
    s += '''mtx <- obj@expr.data\n'''
    s += '''mtx <- t(mtx)         # cell x gene matrix\n'''
    s += '''hc <- hclust(dist(mtx, method = "%s"), method = "%s")\n''' % (dist, hclust)  
    s += '''label <- cutree(tree = hc, k = %d)\n''' % k
    s += '''df <- data.frame(\n'''
    s += '''    barcode = gsub(".", "-", names(label), fixed = TRUE),\n'''
    s += '''    prediction = label - 1\n'''
    s += ''')\n'''
    s += '''write.table(\n'''
    s += '''    df,\n'''
    s += '''    file = "%s",\n''' % out_fn
    s += '''    sep = "\\t",\n'''
    s += '''    row.names = FALSE,\n'''
    s += '''    col.names = TRUE\n'''
    s += ''')\n'''
    s += '''\n'''
    
    script_fn = os.path.join(tmp_dir, "predict_subclones_from_expression.R")
    with open(script_fn, "w") as fp:
        fp.write(s)


    # run the R script.
    if verbose:
        info("run the R script to predict subclones from expression ...")
    exe_cmdline("Rscript %s" % script_fn)
    
    
    # Save to TSV
    info(f"Processed predictions saved to '{out_fn}'.")


    # Print summary
    df = pd.read_csv(out_fn, sep = '\t')
    labels = np.unique(df['prediction'])
    n_cells = len(df)
    info(f"Processed {n_cells} cells after filtering.")
    info("In total %d clones." % len(labels))
    
    return(out_fn)


# UNUSED: Kept for reference. Use predict_subclones_from_expression instead.
def predict_subclones_from_hclust(
    obj_fn,
    out_fn,
    k,
    tmp_dir,
    verbose = False
):
    # Check args.
    assert_e(obj_fn)
    os.makedirs(tmp_dir, exist_ok = True)
    
    
    # Perform subclone identification using existing hierarchical clustering.
    s  = ""
    s += '''# Cut inferCNV hierarchical trees into a specified number of clusters\n'''
    s += '''\n'''
    s += '''obj <- readRDS("%s")\n''' % obj_fn
    s += '''\n'''
    s += '''# Check if clustering data exists\n'''
    s += '''if (is.null(obj@tumor_subclusters$hc)) {\n'''
    s += '''    stop("No hierarchical clustering found. Ensure you ran infercnv::run(analysis_mode='subclusters').")\n'''
    s += '''}\n'''
    s += '''\n'''
    s += '''# Iterate through each group (sample/patient) and cut the tree\n'''
    s += '''all_clusters_list <- lapply(names(obj@tumor_subclusters$hc), function(group_name) {\n'''
    s += '''    \n'''
    s += '''    hc_tree <- obj@tumor_subclusters$hc[[group_name]]\n'''
    s += '''    \n'''
    s += '''    # Perform cuttree\n'''
    s += '''    cuts <- cutree(hc_tree, k = %d)\n''' % k
    s += '''    \n'''
    s += '''    # Create data frame for this group\n'''
    s += '''    data.frame(\n'''
    s += '''        cell_barcode = names(cuts),\n'''
    s += '''        group_id = group_name,\n'''
    s += '''        cluster_id = paste0(group_name, "_c", cuts),\n'''
    s += '''        stringsAsFactors = FALSE\n'''
    s += '''    )\n'''
    s += '''})\n'''
    s += '''\n'''
    s += '''# Combine results from all groups into one data frame\n'''
    s += '''final_df <- do.call(rbind, all_clusters_list)\n'''
    s += '''\n'''
    s += '''# If we have more than k clusters (multiple groups), merge to k using centroid-based hclust\n'''
    s += '''unique_clusters <- unique(final_df$cluster_id)\n'''
    s += '''n_clus <- length(unique_clusters)\n'''
    s += '''if (n_clus > %d) {\n''' % k
    s += '''    # Get expression matrix (cells x genes)\n'''
    s += '''    expr_mtx <- t(obj@expr.data)\n'''
    s += '''    rownames(expr_mtx) <- gsub(".", "-", rownames(expr_mtx), fixed = TRUE)\n'''
    s += '''    expr_mtx <- expr_mtx[rownames(expr_mtx) %in% final_df$cell_barcode, , drop = FALSE]\n'''
    s += '''    # Build centroid matrix (one row per cluster)\n'''
    s += '''    centroid_list <- lapply(unique_clusters, function(cid) {\n'''
    s += '''        cells_c <- final_df$cell_barcode[final_df$cluster_id == cid]\n'''
    s += '''        cells_c <- intersect(cells_c, rownames(expr_mtx))\n'''
    s += '''        if (length(cells_c) == 0) return(rep(NA, ncol(expr_mtx)))\n'''
    s += '''        colMeans(expr_mtx[cells_c, , drop = FALSE])\n'''
    s += '''    })\n'''
    s += '''    centroid_mtx <- do.call(rbind, centroid_list)\n'''
    s += '''    rownames(centroid_mtx) <- unique_clusters\n'''
    s += '''    na_rows <- apply(centroid_mtx, 1, function(x) any(is.na(x)))\n'''
    s += '''    if (any(na_rows)) centroid_mtx <- centroid_mtx[!na_rows, , drop = FALSE]\n'''
    s += '''    if (nrow(centroid_mtx) > 1) {\n'''
    s += '''        hc_centroid <- hclust(dist(centroid_mtx), method = "ward.D2")\n'''
    s += '''        merged <- cutree(hc_centroid, k = %d)\n''' % k
    s += '''        names(merged) <- rownames(centroid_mtx)\n'''
    s += '''        pred_merged <- merged[final_df$cluster_id]\n'''
    s += '''        pred_merged[is.na(pred_merged)] <- 1L\n'''
    s += '''        final_df$prediction <- as.integer(factor(pred_merged)) - 1L\n'''
    s += '''    } else {\n'''
    s += '''        final_df$prediction <- 0L\n'''
    s += '''    }\n'''
    s += '''} else {\n'''
    s += '''    cluster_map <- setNames(0:(n_clus - 1), unique_clusters)\n'''
    s += '''    final_df$prediction <- cluster_map[final_df$cluster_id]\n'''
    s += '''}\n'''
    s += '''\n'''
    s += '''# Format barcode (replace dots with dashes)\n'''
    s += '''final_df$barcode <- gsub(".", "-", final_df$cell_barcode, fixed = TRUE)\n'''
    s += '''\n'''
    s += '''# Select and reorder columns\n'''
    s += '''df <- final_df[, c("barcode", "prediction", "cluster_id")]\n'''
    s += '''\n'''
    s += '''write.table(\n'''
    s += '''    df,\n'''
    s += '''    file = "%s",\n''' % out_fn
    s += '''    sep = "\\t",\n'''
    s += '''    row.names = FALSE,\n'''
    s += '''    col.names = TRUE\n'''
    s += ''')\n'''
    s += '''\n'''
    
    script_fn = os.path.join(tmp_dir, "predict_subclones_from_hclust.R")
    with open(script_fn, "w") as fp:
        fp.write(s)


    # run the R script.
    if verbose:
        info("run the R script to predict subclones from hierarchical clustering ...")
    exe_cmdline("Rscript %s" % script_fn)
    
    
    # Save to TSV
    info(f"Processed predictions saved to '{out_fn}'.")


    # Print summary
    df = pd.read_csv(out_fn, sep = '\t')
    labels = np.unique(df['prediction'])
    n_cells = len(df)
    info(f"Processed {n_cells} cells after filtering.")
    info("In total %d clones." % len(labels))
    
    return(out_fn)



def extract_cna_expression(obj_fn, out_fn, tmp_dir, verbose = False):
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
    def __init__(self, clone_post_fn, run_id = None):
        """Numbat object.
        
        Parameters
        ----------
        clone_post_fn : str
            Path to Numbat TSV file with columns including 'cell' and 
            'clone_opt'.
        run_id : str or None, default None
            Optional run identifier for multiple runs of the same tool.
        """
        super().__init__(tid = "Numbat", run_id = run_id)
        self.clone_post_fn = clone_post_fn
        

    def predict(self, out_fn, verbose = False):
        """Predict subclonal structure from Numbat output.

        Saves a TSV file with columns: `barcode`, `prediction`
        to out_dir/numbat_predictions.tsv.
        """
        return extract_clone_labels(
            clone_post_fn = self.clone_post_fn,
            out_fn = out_fn,
            barcode_col = 'cell',
            clone_col = 'clone_opt',
            delimiter = '\t',
            verbose = verbose
        )



def extract_clone_labels(
    clone_post_fn,
    out_fn,
    barcode_col = 'cell',
    clone_col = 'clone_opt',
    delimiter = '\t',
    verbose = False
):
    # Check args.
    assert_e(clone_post_fn)

    df = pd.read_csv(clone_post_fn, delimiter = delimiter)
    required_cols = [barcode_col, clone_col]
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}")


    # Filter out rows with empty clone labels.
    initial_n_cells = len(df)
    df = df.dropna(subset = [clone_col])
    n_cells = len(df)
    n_removed = initial_n_cells - n_cells
    if n_removed > 0:
        warn("Removed %d cells with empty '%s' values." %  \
            (n_removed, clone_col))
    if n_cells == 0:
        error("No cells remain after removing empty '%s' values." % clone_col)
        raise ValueError


    # Format labels.
    labels_old = np.unique(df[clone_col])
    label_new = np.arange(len(labels_old))
    label_map = {o:n for o, n in zip(labels_old, label_new)}


    # Save to TSV
    result_df = pd.DataFrame({
        'barcode': df[barcode_col].to_numpy(),
        'prediction': df[clone_col].map(label_map),
        'clone_label': df[clone_col]
    })
    result_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Predictions saved to {out_fn}.")


    # Print summary
    df = pd.read_csv(out_fn, sep = '\t')
    labels = np.unique(df['prediction'])
    n_cells = len(df)
    info(f"Processed {n_cells} cells after filtering.")
    info("In total %d clones." % len(labels))

    return(out_fn)



####################################################
#------------------ tools.xclone ------------------#
####################################################

class XClone(Tool):
    def __init__(self, clone_post_fn, run_id = None):
        """XClone object.
        
        Parameters
        ----------
        clone_post_fn : str
            Path to XClone TSV file with columns including 'cell_barcode' and 
            'clone_id_refined'.
        run_id : str or None, default None
            Optional run identifier for multiple runs of the same tool.
        """
        super().__init__(tid = "XClone", run_id = run_id)
        self.clone_post_fn = clone_post_fn
        

    def predict(self, out_fn, verbose = False):
        """Predict subclonal structure from XClone output.

        Saves a TSV file with columns: `barcode`, `prediction`
        to out_dir/xclone_predictions.tsv.
        """
        return extract_xclone_labels(
            clone_post_fn = self.clone_post_fn,
            out_fn = out_fn,
            barcode_col = 'cell_barcode',
            clone_col = 'clone_id_refined',
            delimiter = '\t',
            verbose = verbose
        )



def extract_xclone_labels(
    clone_post_fn,
    out_fn,
    barcode_col = 'cell_barcode',
    clone_col = 'clone_id_refined',
    delimiter = '\t',
    verbose = False
):
    # Check args.
    assert_e(clone_post_fn)

    df = pd.read_csv(clone_post_fn, delimiter = delimiter)
    
    # Try to find barcode column if default doesn't exist
    if barcode_col not in df.columns:
        # Try common alternatives
        for alt_col in ['cell_barcode', 'cell', 'Cell', 'CELL', 'barcode', 'Barcode', 'BARCODE']:
            if alt_col in df.columns:
                barcode_col = alt_col
                if verbose:
                    info(f"Using '{barcode_col}' as barcode column.")
                break
        else:
            raise ValueError(f"TSV must contain a barcode column. Found columns: {list(df.columns)}")
    
    required_cols = [barcode_col, clone_col]
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"TSV must contain columns: {required_cols}. Found columns: {list(df.columns)}")


    # Filter out rows with empty clone labels.
    initial_n_cells = len(df)
    df = df.dropna(subset = [clone_col])
    n_cells = len(df)
    n_removed = initial_n_cells - n_cells
    if n_removed > 0:
        warn("Removed %d cells with empty '%s' values." %  \
            (n_removed, clone_col))
    if n_cells == 0:
        error("No cells remain after removing empty '%s' values." % clone_col)
        raise ValueError


    # Format labels.
    labels_old = np.unique(df[clone_col])
    label_new = np.arange(len(labels_old))
    label_map = {o:n for o, n in zip(labels_old, label_new)}


    # Save to TSV
    result_df = pd.DataFrame({
        'barcode': df[barcode_col].to_numpy(),
        'prediction': df[clone_col].map(label_map),
        'clone_label': df[clone_col]
    })
    result_df.to_csv(out_fn, sep = '\t', index = False)
    info(f"Predictions saved to {out_fn}.")


    # Print summary
    df = pd.read_csv(out_fn, sep = '\t')
    labels = np.unique(df['prediction'])
    n_cells = len(df)
    info(f"Processed {n_cells} cells after filtering.")
    info("In total %d clones." % len(labels))

    return(out_fn)




##################################################
#------------------ utils.base ------------------#
##################################################

def assert_e(path):
    """Assert file or folder exists, mimicking shell "test -e"."""
    assert path is not None and os.path.exists(path)

    
def assert_n(x):
    """Assert `x` has content, mimicking shell "test -n"."""
    assert x is not None and len(x) > 0
    
    
    
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
