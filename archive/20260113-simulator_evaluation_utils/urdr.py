# urdr.py


import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import statsmodels.api as sm

from udf import long2wide
from umath import get_log2FC
from umatrix import sparse2array
from uplot import ax_format_wrapper, savefig_wrapper

from uplot import box_plot as u_box_plot
from uplot import stat_LR as u_stat_LR
from uplot import scatter_plot_LR as u_scatter_plot_LR
from uplot import mulax_LR as u_mulax_LR
from uplot import stack_plot_perc as u_stack_plot_perc
from uplot import violin_plot as u_violin_plot
from usmooth import stat_smooth as u_stat_smooth
from usmooth import smooth_plot as u_smooth_plot



def box_plot(
    x, y, ylabel, data = None, log_scale = False, 
    colors = None, showfliers = True, rotation = 45, fill = False,
    ylims = None,
    ax_kws = None,
    savefig_kws = None
):
    """Box plot for cell-wise or gene-wise metrics."""
    return u_box_plot(
        x = x, y = y, ylabel = ylabel, data = data, log_scale = log_scale,
        colors = colors, showfliers = showfliers, rotation = rotation, fill = fill,
        ylims = ylims,
        ax_kws = ax_kws,
        savefig_kws = savefig_kws
    )



def stat_log2fc(
    data,
    metric = "mean", group = "X_name", gvars = None
):
    """log2FC.
    gvars (None or DataFrame).
    """
    df = long2wide(data, columns = group, values = metric)

    if gvars is None:
        groups = data[group].unique()
        cx, cy = [], []
        n = len(groups)
        for i in range(n):
            for j in range(i):
                cx.append(groups[j])
                cy.append(groups[i])
        gvars = pd.DataFrame(data = dict(x = cx, y = cy))
    else:
        assert 'x' in gvars.columns
        assert 'y' in gvars.columns

    res = None
    for i in range(gvars.shape[0]):
        cx = gvars['x'].loc[i]
        cy = gvars['y'].loc[i]
        x = df[cx].to_numpy()
        y = df[cy].to_numpy()
        d = pd.DataFrame(data = {
                'x': cx,
                'y': cy,
                metric + '_x': x,
                metric + '_y': y,
                'log2FC': get_log2FC(x, y)
        })
        if i == 0:
            res = d
        else:
            res = pd.concat([res, d], ignore_index = True)
    return(res)



@ax_format_wrapper
def __log2fc_vs_mean(
    ax,
    x, y, data,
    xlabel = None, ylabel = None,
    log_scale_xaxis = False, 
    cut = None, cut_markers = None,
    point_size = None, point_color = None,
    fc_line = True, fc_kws = None,
    ref_line = True, ref_kws = None,
    ylims = None,
    ax_kws = None
):
    """Scatter plot for log2FC of mean (y-axis) vs. mean (x-axis).
    refer to the DESeq2 paper.
    x, y (str): column names.
    data (DataFrame): wide-format.
    log_scale_xaxis (bool): whether to transform x axis with log10(i+1).
    cut (None or tuple of 2 float): low and high cutoff of y axis.
    cut_markers (None or tuple of 2 str): scatter markers for low and up 
        cutoffs.
    point_size (None or float): scatter point size.
    fc_line (bool): whether to plot the FC lines (`y = -1` and `y = 1`).
    ref_line (bool): whether to plot the reference line (`y = 0`).
    """
    cx, cy = x, y
    vx = data[cx].to_numpy()
    vy = data[cy].to_numpy()
    if log_scale_xaxis:
        vx = np.log10(vx + 1)

    if cut is None:
        ax.plot(vx, vy, marker = '.', markersize = point_size, 
                linestyle = 'none', color = point_color)
    else:
        assert len(cut) == 2
        assert len(cut_markers) == 2
        low, high = cut
        
        x1 = vx[vy > high]
        y1 = np.repeat(high, len(x1))
        x2 = vx[np.logical_and(vy >= low, vy <= high)]
        y2 = vy[np.logical_and(vy >= low, vy <= high)]
        x3 = vx[vy < low]
        y3 = np.repeat(low, len(x3))
    
        ax.plot(x1, y1, marker = cut_markers[1], markersize = point_size, 
                linestyle = 'none', color = point_color)
        ax.plot(x2, y2, marker = '.', markersize = point_size, 
                linestyle = 'none', color = point_color)
        ax.plot(x3, y3, marker = cut_markers[0], markersize = point_size, 
                linestyle = 'none', color = point_color)
        
    if fc_line:
        if fc_kws is None:
            fc_kws = dict()
        ax.axhline(
            y = fc_kws.get('y_up', 1),
            linestyle = fc_kws.get('linestyle', '--'),
            linewidth = fc_kws.get('linewidth', 1.0),
            color = fc_kws.get('color', 'black')
        )
        ax.axhline(
            y = fc_kws.get('y_low', -1),
            linestyle = fc_kws.get('linestyle', '--'),
            linewidth = fc_kws.get('linewidth', 1.0),
            color = fc_kws.get('color', 'black')
        )
    
    if ref_line:
        if ref_kws is None:
            ref_kws = dict()
        ax.axhline(
            y = ref_kws.get('y', 0), 
            linestyle = ref_kws.get('linestyle', '--'), 
            linewidth = ref_kws.get('linewidth', 1.0),
            color = ref_kws.get('color', 'gray')
        )
        
    if ylims is not None:
        ax.set_ylim(*ylims)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return(ax)



@savefig_wrapper
def scatter_plot_log2fc_vs_mean(
    x, y, data,
    title = None,
    xlabel = None, ylabel = None,
    log_scale_xaxis = False, 
    cut = None, cut_markers = None,
    point_size = None, point_color = None,
    fc_line = True, fc_kws = None,
    ref_line = True, ref_kws = None,
    ylims = None,
    ax_kws = None,
    savefig_kws = None
):
    """Scatter plot for log2FC of mean (y-axis) vs. mean (x-axis).
    refer to the DESeq2 paper.
    x, y (str): column names.
    data (DataFrame): wide-format.
    log_scale_xaxis (bool): whether to transform x axis with log10(i+1).
    cut (None or tuple of 2 float): low and high cutoff of y axis.
    cut_markers (None or tuple of 2 str): scatter markers for low and up 
        cutoffs.
    point_size (None or float): scatter point size.
    ref_line_color (None or color): color for reference line (y = 0).
    """
    ax = plt.gca()
    __log2fc_vs_mean(
        ax = ax,
        x = x, y = y, data = data,
        xlabel = xlabel, ylabel = ylabel,
        log_scale_xaxis = log_scale_xaxis, 
        cut = cut, cut_markers = cut_markers,
        point_size = point_size, point_color = point_color,
        fc_line = fc_line, fc_kws = fc_kws,
        ref_line = ref_line, ref_kws = ref_kws,
        ylims = ylims,
        ax_kws = ax_kws
    )
    ax.set_title(title)
    return(ax)



@savefig_wrapper
def mulax_log2fc_vs_mean(
    data,
    x, y, vx, vy,
    xlabel = None, ylabel = None,
    log_scale_xaxis = False,
    cut = None, cut_markers = None,
    point_size = None, point_color = None,
    fc_line = True, fc_kws = None,
    ref_line = True, ref_kws = None,
    ylims = None,
    xlabel_pos = None, ylabel_pos = None,
    title_font_size = 10,
    nrows = 1, ncols = 1,
    ax_kws = None,
    savefig_kws = None
):
    """Multiple axes for log2FC vs. gene mean.
    gvars (None or DataFrame).
    """
    gvars = data[[x, y]].drop_duplicates()
    gvars = gvars.rename(columns = {x: 'x', y: 'y'})
    gvars["xy_group"] = gvars['y'] + ' ~ ' + gvars['x']

    g = sns.FacetGrid(gvars, col = 'xy_group', col_wrap = ncols)
    
    k = 0
    for i in range(gvars.shape[0]):
        ax = g.axes[k]
        cx = gvars['x'].iloc[i]
        cy = gvars['y'].iloc[i]
        df = data.loc[(data[x] == cx) & (data[y] == cy)].copy()
        __log2fc_vs_mean(
            ax = ax,
            x = vx, y = vy, data = df,
            xlabel = None, ylabel = None,
            log_scale_xaxis = log_scale_xaxis,
            cut = cut, cut_markers = cut_markers,
            point_size = point_size, point_color = point_color,
            fc_line = fc_line, fc_kws = fc_kws,
            ref_line = ref_line, ref_kws = ref_kws,
            ylims = ylims,
            ax_kws = ax_kws
        )
        ax.set_title(cy + ' ~ ' + cx, fontsize = title_font_size)
        k += 1
        
    g.set_xlabels(label = xlabel, clear_inner = True)
    if xlabel_pos is not None:
        for ax in g.axes:
            if ax.get_xlabel():
                ax.xaxis.set_label_coords(xlabel_pos[0], xlabel_pos[1])
    
    g.set_ylabels(label = ylabel, clear_inner = True)
    if ylabel_pos is not None:
        for ax in g.axes:
            if ax.get_ylabel():
                ax.yaxis.set_label_coords(ylabel_pos[0], ylabel_pos[1])

    for ax in g.axes:
        for s in ["top", "right"]:
            ax.spines[s].set_visible(True)
    return(g)



def stat_LR(
    data,
    metric = "mean", group = "X_name", gvars = None, 
    log_scale = False
):
    return u_stat_LR(
        data = data,
        metric = metric, group = group, gvars = gvars, 
        log_scale = log_scale
    )



def scatter_plot_LR(
    x, y, data, 
    lr_kws,
    pos_r, pos_fit,
    point_size,
    title = None,
    xlabel = None, ylabel = None,
    log_scale = False,
    fit_line = True, fit_kws = None,
    ref_line = False, ref_kws = None,
    xlims_line = None,
    ylims = None,
    ax = None,
    ax_kws = None,
    savefig_kws = None
):
    """Scatter plot for linear regression."""
    return u_scatter_plot_LR(
        x = x, y = y, data = data,
        lr_kws = lr_kws,
        pos_r = pos_r, pos_fit = pos_fit,
        point_size = point_size, 
        title = title,
        xlabel = xlabel, ylabel = ylabel,
        log_scale = log_scale,
        fit_line = fit_line, fit_kws = fit_kws,
        ref_line = ref_line, ref_kws = ref_kws,
        xlims_line = xlims_line,
        ylims = ylims,
        ax = ax, ax_kws = ax_kws,
        savefig_kws = savefig_kws
    )
    


def mulax_LR(
    data, st_lr, 
    pos_r, pos_fit, 
    columns = None,
    xlabel = None, ylabel = None,
    log_scale = False,
    point_size = 25,
    fit_line = True, fit_kws = None,
    ref_line = False, ref_kws = None,
    ylims = None,
    xlabel_pos = None, ylabel_pos = None,
    title_font_size = 10,
    nrows = 1, ncols = 1,
    ax_kws = None,
    savefig_kws = None
):
    """Multiple axes for linear regressions."""
    return u_mulax_LR(
        data = data, st_lr = st_lr,
        pos_r = pos_r, pos_fit = pos_fit, 
        columns = columns, 
        xlabel = xlabel, ylabel = ylabel,
        log_scale = log_scale,
        point_size = point_size,
        fit_line = fit_line, fit_kws = fit_kws,
        ref_line = ref_line, ref_kws = ref_kws,
        ylims = ylims,
        xlabel_pos = xlabel_pos, ylabel_pos = ylabel_pos,
        title_font_size = title_font_size,
        nrows = nrows, ncols = ncols,
        ax_kws = ax_kws,
        savefig_kws = savefig_kws
    )



def stat_smooth(
    x, y, data, 
    group = "X_name", 
    log_scale = (False, False),
    smooth_kws = None
):
    return u_stat_smooth(
        x = x, y = y, data = data, 
        group = group, 
        log_scale = log_scale,
        smooth_kws = smooth_kws
    )
            
            

@savefig_wrapper
@ax_format_wrapper
def smooth_plot(
    x, y, data, 
    group = "X_name", 
    xlabel = None, ylabel = None,
    legend = False, legend_pos = None,
    line_width = 1,
    ylims = None,
    ax_kws = None,
    savefig_kws = None
):
    return u_smooth_plot(
        x = x, y = y, data = data, 
        group = group, 
        xlabel = xlabel, ylabel = ylabel,
        legend = legend, legend_pos = legend_pos,
        line_width = line_width,
        ylims = ylims,
        ax_kws = ax_kws,
        savefig_kws = savefig_kws
    )



def __calc_gene_percent(df, metric, cutoffs, labels):
    """Calculate the percentage of genes in specific range of zero 
    proportions."""
    values = []
    cumsum = 0
    for i, c in enumerate(cutoffs):
        v = None
        if i < len(cutoffs) - 1:
            v = (df[metric] < c).mean() * 100
        else:
            v = (df[metric] <= c).mean() * 100
        v -= cumsum
        cumsum += v
        values.append(v)
    res = pd.DataFrame(
        data = dict(cutoff = cutoffs, label = labels, count = values))
    return(res)



def stat_stack_zeroprop(data, metric = "zero_prop", group = "X_name", 
                        cutoffs = None, labels = None):
    """Stack bar plot for percentage of genes within various ranges of 
    gene zero proportions.
    data (DataFrame).
    metric (str): column name in `data` used as values of groups Xs.
    group (str): column name in `data`. That column contains group names.
    cutoffs (None or list of float): a list of cutoff values.
    labels (None or list of str): labels for `cutoffs`.
    """
    if cutoffs is None:
        assert labels is None
        cutoffs = [.8, .9, .95, .99, 1]     #cutoffs = [.9, .95, .99, 1, 1.01]
    else:
        assert labels is not None
    if labels is None:
        labels = ["[0, 0.8)", "[0.8, 0.9)", "[0.9, 0.95)", "[0.95, 0.99)", "[0.99, 1]"]
    assert len(cutoffs) == len(labels)
    res = data.groupby(group).apply(
        __calc_gene_percent, 
        metric = metric, cutoffs = cutoffs, labels = labels
    )
    res = res.loc[data[group].unique()].reset_index()
    res = res.rename(columns = dict(count = 'percent'))
    return(res)



def stack_plot_perc(
    x, y, data, hue, 
    ylabel, legend_title, 
    font_size, bbox_to_anchor, colors,
    rotation = 45,
    min_v = 0,
    ax_kws = None,
    savefig_kws = None
):
    """Stacked bar plot for percentage of groups."""
    return u_stack_plot_perc(
        x = x, y = y, data = data, hue = hue, 
        ylabel = ylabel, legend_title = legend_title, 
        font_size = font_size, bbox_to_anchor = bbox_to_anchor, 
        colors = colors,
        rotation = rotation,
        min_v = min_v,
        ax_kws = ax_kws,
        savefig_kws = savefig_kws
    )



def violin_plot(
    x, y, ylabel, data = None, log_scale = False, 
    colors = None, rotation = 45, inner = None, cut = 2, fill = False,
    ylims = None,
    ax_kws = None,
    savefig_kws = None
):
    """Violin plot for cell-wise or gene-wise metrics."""
    return u_violin_plot(
        x = x, y = y, ylabel = ylabel, data = data, 
        log_scale = log_scale, 
        colors = colors, rotation = rotation, inner = inner, 
        cut = cut, fill = fill,
        ylims = ylims,
        ax_kws = ax_kws,
        savefig_kws = savefig_kws
    )



### Metrics ###


### Cell-wise metrics


def __get_metrics(mtype, X, metrics = None, out_fmt = "df"):
    """Wrapper function for metrics calculation.

    Parameters
    ----------
    mtype : {"cw", "gw"}
        Metric type. One of "cw" (cell-wise) or "gw" (gene-wise).
    X : matrix-like
        The *cell x feature* matrix.
    metrics : list of str or None, default None
        A list of metrics to be calculated, each of which should be among
        "lib_size" (library size), "zero_prop" (zero proportion) if cell-wise;
        or "mean" (mean), "var" (variance), "cv" (coefficient of variation),
        "zero_prop" (zero proportion) if gene-wise.
        If None, all available metrics will be used.
    out_fmt : {"df", "dict", "list"}
        Format of the returned result. 
        One of "df" (pandas.DataFrame), "dict" (dict), and "list" (list).

    Returns
    -------
    object
        An object in the format specified by `out_fmt`.
    """
    assert mtype in ("cw", "gw")

    all_metrics = None
    if mtype == "cw":
        all_metrics = CELLWISE_METRICS
    else:
        all_metrics = GENEWISE_METRICS

    all_out_fmt = ["df", "dict", "list"]

    X = sparse2array(X)

    if metrics is None:
        metrics = all_metrics
    for m in metrics:
        if m not in all_metrics:
            error("invalid metric '%s'." % m)
            raise ValueError
        
    if out_fmt not in all_out_fmt:
        error("invalid output format '%s'." % out_fmt)
        raise ValueError
    
    res = []
    if mtype == "cw":
        for m in metrics:
            if m == "lib_size":
                res.append(get_cw_lib_size(X))
            elif m == "zero_prop":
                res.append(get_cw_zero_prop(X))
            else:
                error("invalid metric '%s'." % m)
                raise ValueError
    else:
        for m in metrics:
            if m == "mean":
                res.append(get_gw_mean(X))
            elif m == "var":
                res.append(get_gw_var(X))
            elif m == "cv":
                res.append(get_gw_cv(X))
            elif m == "zero_prop":
                res.append(get_gw_zero_prop(X))
            else:
                error("invalid metric '%s'." % m)
                raise ValueError
        
    if out_fmt == "df":
        res = pd.DataFrame(data = {m:v for m, v in zip(metrics, res)})
    elif out_fmt == "dict":
        res = {m:v for m, v in zip(metrics, res)}
    elif out_fmt == "list":
        pass
    else:
        error("invalid output format '%s'." % out_fmt)
        raise ValueError
    
    return(res)


def __get_metrics_group(mtype, X_lst, id_lst, X_names = None, metrics = None):
    """Wrapper function for metrics calculation in a group of matrices.

    Parameters
    ----------
    mtype : {"cw", "gw"}
        Metric type. One of "cw" (cell-wise) or "gw" (gene-wise).
    X_lst : list of matrix-like
        A list of *cell x feature* matrices (matrix-like objects).
    X_names : list of str or None, default None
        A list of group names (str).
        Its length and order should match `X_lst`.
        If `None`, the default ["X0", "X1", ..., "Xn"] will be used.
    metrics : list of str or None, default None
        A list of metrics to be calculated, each of which should be among
        "lib_size" (library size), "zero_prop" (zero proportion) if cell-wise;
        or "mean" (mean), "var" (variance), "cv" (coefficient of variation),
        "zero_prop" (zero proportion) if gene-wise.
        If None, all available metrics will be used.

    Returns
    -------
    pandas.DataFrame
        A `pandas.DataFrame` object containing calculated metrics, whose
        first several column names are the `metrics` and the last column is 
        "X_name" storing the group names.
    """
    assert mtype in ("cw", "gw")

    if X_names is None:
        X_names = ["X" + str(i) for i in range(len(X_lst))]
    if len(X_lst) != len(X_names):
        error("length of 'X_lst' and 'X_names' should be the same!")
        raise ValueError
    
    result = None
    for idx, (X, ids, name) in enumerate(zip(X_lst, id_lst, X_names)):
        res = __get_metrics(mtype, X, metrics = metrics, out_fmt = "df")
        res["X_name"] = name
        if mtype == 'cw':
            res['cell'] = ids
        else:
            res['feature'] = ids
        if idx == 0:
            result = res
        else:
            result = pd.concat([result, res], ignore_index = True)

    return(result)


def get_cw_metrics_group(X_lst, id_lst, X_names = None, metrics = None):
    return(__get_metrics_group(
        mtype = "cw",
        X_lst = X_lst,
        id_lst = id_lst,
        X_names = X_names,
        metrics = metrics
    ))


def get_cw_metrics(X, metrics = None, out_fmt = "df"):
    return(__get_metrics(
        mtype = "cw",
        X = X,
        metrics = metrics,
        out_fmt = out_fmt
    ))


def get_cw_lib_size(X):
    return np.sum(X, axis = 1)

def get_cw_zero_prop(X):
    return np.mean(X <= 0.0, axis = 1)


CELLWISE_METRICS = ["lib_size", "zero_prop"]



### Gene-wise metrics

def get_gw_metrics_group(X_lst, id_lst, X_names = None, metrics = None):
    return(__get_metrics_group(
        mtype = "gw",
        X_lst = X_lst,
        id_lst = id_lst,
        X_names = X_names,
        metrics = metrics
    ))


def get_gw_metrics(X, metrics = None, out_fmt = "df"):
    return(__get_metrics(
        mtype = "gw",
        X = X,
        metrics = metrics,
        out_fmt = out_fmt
    ))
  

def get_gw_cv(X):
    """Get gene-wise coefficient of variation (CV)."""
    return np.std(X, axis = 0) / (np.mean(X, axis = 0) + 1e-8)

def get_gw_mean(X):
    """Get gene-wise mean."""
    return np.mean(X, axis = 0)

def get_gw_var(X):
    """Get gene-wise variance."""
    return np.var(X, axis = 0)

def get_gw_zero_prop(X):
    """Get gene-wise zero proportion."""
    return np.mean(X <= 0.0, axis = 0)


GENEWISE_METRICS = ["mean", "var", "cv", "zero_prop"]
