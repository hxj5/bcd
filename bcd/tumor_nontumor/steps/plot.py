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
        tool_col = 'tool',
        tool_list = tool_list,
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
    if fig_nrow is None:
        fig_nrow = np.ceil(n / fig_ncol)
    if fig_width is None:
        fig_width = 4 * ncol
    if fig_height is None:
        fig_height = 4 * fig_nrow
    fig, axes = plt.subplots(
        fig_nrow, fig_ncol, 
        figsize = (fig_width, fig_height)
    )
    if n == 1:
        axes = [axes]

    df = pd.read_csv(truth_fn, sep = '\t')
    truth_lables = df['annotation'].to_numpy()
    for ax, tool, tool_fn in zip(axes, tool_list, tool_fn_list):
        tid = tool.tid
        if verbose:
            info("process %s ..." % tid)

        df = pd.read_csv(tool_fn, sep = '\t')
        tool_lables = df['prediction'].to_numpy()
        cm = confusion_matrix(
            truth_lables, tool_lables, 
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
    fig_width = 10,
    fig_height = 5,
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
        t = (df['Prediction'] == 'tumor').sum()
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
    tool_col = 'tool',
    tool_list = None,
    metric_list = None,
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
    assert tool_col in df.columns
    if tool_list is None:
        tool_list = df[tool_col].to_numpy()
    else:
        for tool in tool_list:
            assert tool.tid in df[tool_col]
    if metric_list is None:
        s = df.columns[df.columns != tool_tol]
        if len(s) <= 0:
            raise ValueError("no metric columns in '%s'." % metric_fn)
        metric_list = s.to_numpy()
    else:
        for metric in metric_list:
            assert metric in df.columns


    # plot bar for each tool.
    if verbose:
        info("plot bar for each tool's metrics ...")
    
    n = len(tool_list)
    if fig_nrow is None:
        fig_nrow = np.ceil(n / fig_ncol)
    if fig_width is None:
        fig_width = 4 * ncol
    if fig_height is None:
        fig_height = 5 * fig_nrow
    fig, axes = plt.subplots(
        fig_nrow, fig_ncol,
        figsize = (fig_width, fig_height),
        sharey = True
    )
    if n == 1:
        axes = [axes]    

    bar_width = 0.6
    pos       = np.arange(len(tool_list))

    for ax, metric in zip(axes, metric_list):
        values = [df.loc[df[tool_col] == tool.tid, metric] \
                  for tool in tool_list]
        bars = ax.bar(pos, values, bar_width, 
                      color = '#1f77b4', edgecolor = 'black')
        ax.set_title(metric, fontsize = 12)
        ax.set_xticks(pos)
        ax.set_xticklabels([tool.tid for tool in tool_list], 
                           rotation = 45, ha = 'right')
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
    tool_col = 'tool',
    tool_list = None,
    metric_list = None,
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
    assert tool_col in df.columns
    if tool_list is None:
        tool_list = df[tool_col].to_numpy()
    else:
        for tool in tool_list:
            assert tool.tid in df[tool_col]
    if metric_list is None:
        s = df.columns[df.columns != tool_tol]
        if len(s) <= 0:
            raise ValueError("no metric columns in '%s'." % metric_fn)
        metric_list = s.to_numpy()
    else:
        for metric in metric_list:
            assert metric in df.columns


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
    for tool in tool_list:
        tid = tool.tid
        if verbose:
            info("process %s ..." % tid)

        values = [df.loc[df[tool_col] == tid, m] for m in metric_list]
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
