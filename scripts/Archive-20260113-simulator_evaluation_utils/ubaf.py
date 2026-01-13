# ubaf.py


import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import statsmodels.api as sm

from udf import long2wide
from umatrix import sparse2array
from uplot import ax_format_wrapper, savefig_wrapper

from uplot import box_plot as u_box_plot
from uplot import stat_LR as u_stat_LR
from uplot import scatter_plot_LR as u_scatter_plot_LR
from uplot import mulax_LR as u_mulax_LR
from uplot import stack_plot_perc as u_stack_plot_perc
from uplot import violin_plot as u_violin_plot
from urdr import stat_stack_zeroprop as u_stat_stack_zeroprop
from usmooth import smooth_plot as u_smooth_plot



@savefig_wrapper
@ax_format_wrapper
def box_plot(
    x, y, ylabel, data = None, log_scale = False, 
    colors = None, showfliers = True, rotation = 45, fill = False,
    ylims = None,
    add_scatter = False,
    scatter_kws = None,
    ax_kws = None,
    savefig_kws = None
):
    """Box plot for cell-wise or gene-wise metrics."""
    ax = u_box_plot(
        x = x, y = y, ylabel = ylabel, data = data, log_scale = log_scale,
        colors = colors, showfliers = showfliers, rotation = rotation, fill = fill,
        ylims = ylims,
        ax_kws = None,
        savefig_kws = None
    )
    if add_scatter:
        if data is not None:
            x = data[x]
            y = data[y]
        if log_scale:
            y = np.log10(y + 1)
        if scatter_kws is None:
            sns.stripplot(data = data, x = x, y = y, palette = colors, 
                      hue = x, legend = False, ax = ax)
        else:
            if 'color' in scatter_kws:
                sns.stripplot(data = data, x = x, y = y, 
                      legend = False, ax = ax, **scatter_kws)
            else:
                sns.stripplot(data = data, x = x, y = y, palette = colors,
                      hue = x, legend = False, ax = ax, **scatter_kws)                
    ax.set_xlabel(None)
    return(ax)



@ax_format_wrapper
def __baf_vs_mean(
    ax,
    x, y, data,
    xlabel = None, ylabel = None,
    log_scale_xaxis = False,
    point_size = None, point_color = None,
    ref_line = True, ref_kws = None,
    ylims = None,
    ax_kws = None
):
    """Scatter plot for BAF (y-axis) vs. mean (x-axis).
    x, y (str): column names.
    data (DataFrame): wide-format.
    log_scale_xaxis (bool): whether to transform x axis with log10(i+1).
    point_size (None or float): scatter point size.
    point_color (None or color) scatter point color.
    ref_line (bool): whether to plot the reference line (`y = 0.5`).
    """
    cx, cy = x, y
    x = data[cx].to_numpy()
    y = data[cy].to_numpy()
    if log_scale_xaxis:
        x = np.log10(x + 1)

    ax.plot(x, y, marker = '.', markersize = point_size, 
            linestyle = 'none', color = point_color)
    
    if ref_line:
        if ref_kws is None:
            ref_kws = dict()
        ax.axhline(
            y = ref_kws.get('y', 0.5), 
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
def scatter_plot_baf_vs_mean(
    x, y, data,
    title = None,
    xlabel = None, ylabel = None,
    log_scale_xaxis = False,
    point_size = None, point_color = None,
    ref_line = True, ref_kws = None,
    ylims = None,
    ax_kws = None,
    savefig_kws = None
):
    """Scatter plot for BAF (y-axis) vs. mean (x-axis).
    x, y (str): column names.
    data (DataFrame): wide-format.
    log_scale_xaxis (bool): whether to transform x axis with log10(i+1).
    point_size (None or float): scatter point size.
    point_color (None or color): scatter point color.
    ref_line (bool): whether to plot the reference line (`y = 0.5`).
    """
    ax = plt.gca()
    ax = __baf_vs_mean(
        ax = ax,
        x = x, y = y, data = data,
        xlabel = xlabel, ylabel = ylabel,
        log_scale_xaxis = log_scale_xaxis,
        point_size = point_size, point_color = point_color,
        ref_line = ref_line, ref_kws = ref_kws,
        ylims = ylims,
        ax_kws = ax_kws
    )
    ax.set_title(title)
    return(ax)



@savefig_wrapper
def mulax_baf_vs_mean(
    data,
    x, y, group = "X_name",
    xlabel = None, ylabel = None,
    log_scale_xaxis = False,
    point_size = None, point_color = None,
    ref_line = True, ref_kws = None,
    ylims = None,
    xlabel_pos = None, ylabel_pos = None,
    title_font_size = 10,
    nrows = 1, ncols = 1,
    ax_kws = None,
    savefig_kws = None
):
    """Multiple axes for BAF vs. gene DP mean.
    """
    groups = data[group].unique()
    gvars = pd.DataFrame(data = dict(xy_group = groups))
    g = sns.FacetGrid(gvars, col = 'xy_group', col_wrap = ncols)
    
    k = 0
    for i in range(groups.shape[0]):
        ax = g.axes[k]
        gname = groups[i]
        df = data.loc[data[group] == gname].copy()
        __baf_vs_mean(
            ax = ax,
            x = x, y = y, data = df,
            xlabel = xlabel, ylabel = ylabel,
            log_scale_xaxis = log_scale_xaxis,
            point_size = point_size, point_color = point_color,
            ref_line = ref_line, ref_kws = ref_kws,
            ylims = ylims,
            ax_kws = ax_kws
        )
        ax.set_title(gname, fontsize = title_font_size)
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



@ax_format_wrapper
def __pair_baf_vs_mean(
    ax,
    x, y, data, group, gvars = None,
    xlabel = None, ylabel = None,
    log_scale_xaxis = False,
    point_size = None, point_colors = None,
    ref_line = True, ref_kws = None,
    ylims = None,
    legend_kws = None,
    ax_kws = None
):
    """Scatter plot for pair-wise BAF (y-axis) vs. mean (x-axis).
    x, y (str): column names. Metrics.
    data (DataFrame): long-format.
    group (str): column names. Groups.
    gvars (list of str): a list of groups to plot.
    log_scale_xaxis (bool): whether to transform x axis with log10(i+1).
    point_size (None or float): scatter point size.
    point_colors (list of colors): its order matches `gvars`.
    ref_line (bool): whether to plot the reference line (`y = 0.5`).
    """
    assert len(gvars) == len(point_colors)
    data = data.copy()
    cx, cy = x, y
    if log_scale_xaxis:
        data[cx] = np.log10(data[cx] + 1)

    for g, c in zip(gvars, point_colors):
        df = data.loc[data[group] == g].copy()
        x = df[cx].to_numpy()
        y = df[cy].to_numpy()
        ax.plot(x, y, marker = '.', markersize = point_size, 
                linestyle = 'none', color = c, label = g)

    if ref_line:
        if ref_kws is None:
            ref_kws = dict()
        ax.axhline(
            y = ref_kws.get('y', 0.5), 
            linestyle = ref_kws.get('linestyle', '--'), 
            linewidth = ref_kws.get('linewidth', 1.0),
            color = ref_kws.get('color', 'gray')
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if ylims is not None:
        ax.set_ylim(*ylims)
    if legend_kws is not None:
        ax.legend(**legend_kws)
    return(ax)



@savefig_wrapper
def scatter_plot_pair_baf_vs_mean(
    x, y, data, group, gvars = None,
    title = None,
    xlabel = None, ylabel = None,
    log_scale_xaxis = False,
    point_size = None, point_colors = None,
    ref_line = True, ref_kws = None,
    ylims = None,
    legend_kws = None,
    ax_kws = None,
    savefig_kws = None
):
    """Scatter plot for pair-wise BAF (y-axis) vs. mean (x-axis).
    """
    ax = plt.gca()
    ax = __pair_baf_vs_mean(
        ax = ax,
        x = x, y = y, data = data,
        group = group, gvars = gvars,
        xlabel = xlabel, ylabel = ylabel,
        log_scale_xaxis = log_scale_xaxis,
        point_size = point_size, point_colors = point_colors,
        ref_line = ref_line, ref_kws = ref_kws,
        ylims = ylims,
        legend_kws = legend_kws,
        ax_kws = ax_kws
    )
    return(ax)



@savefig_wrapper
def mulax_pair_baf_vs_mean(
    x, y, data, group = "X_name", gvars = None,
    xlabel = None, ylabel = None,
    log_scale_xaxis = False,
    point_size = None, point_colors = None,
    ref_line = True, ref_kws = None,
    ylims = None,
    legend_kws = None,
    xlabel_pos = None, ylabel_pos = None,
    title_font_size = 10,
    nrows = 1, ncols = 1,
    ax_kws = None,
    savefig_kws = None
):
    """Multiple axes for BAF vs. gene DP mean.
    gvars (None or DataFrame).
    """
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
    gvars["xy_group"] = gvars['y'] + ' ~ ' + gvars['x']

    g = sns.FacetGrid(gvars, col = 'xy_group', col_wrap = ncols)
    
    k = 0
    for i in range(gvars.shape[0]):
        ax = g.axes[k]
        __pair_baf_vs_mean(
            ax = ax,
            x = x, y = y, data = data,
            group = group,
            gvars = [gvars['x'].iloc[i], gvars['y'].iloc[i]],
            xlabel = xlabel, ylabel = ylabel,
            log_scale_xaxis = log_scale_xaxis,
            point_size = point_size, point_colors = point_colors,
            ref_line = ref_line, ref_kws = ref_kws,
            ylims = ylims,
            legend_kws = legend_kws,
            ax_kws = ax_kws
        )
        #ax.set_title(gvars['y'].loc[i] + ' ~ ' + gvars['x'].loc[i], 
        #             fontsize = title_font_size)
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
            
    g.set_titles(col_template = "")
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
    return u_stat_stack_zeroprop(
        data = data, metric = metric, group = group, 
        cutoffs = cutoffs, labels = labels
    )



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
    ax = u_stack_plot_perc(
        x = x, y = y, data = data, hue = hue, 
        ylabel = ylabel, legend_title = legend_title, 
        font_size = font_size, bbox_to_anchor = bbox_to_anchor, 
        colors = colors,
        rotation = rotation,
        min_v = min_v,
        ax_kws = ax_kws,
        savefig_kws = savefig_kws
    )
    return(ax)



@savefig_wrapper
@ax_format_wrapper
def violin_plot(
    x, y, ylabel, data = None, log_scale = False, 
    colors = None, rotation = 45, inner = None, cut = 2, fill = False,
    fontsize = None,
    ylims = None,
    ax_kws = None,
    savefig_kws = None
):
    """Violin plot for cell-wise or gene-wise metrics."""
    ax = u_violin_plot(
        x = x, y = y, ylabel = ylabel, data = data, 
        log_scale = log_scale, 
        colors = colors, rotation = rotation, inner = inner, 
        cut = cut, fill = fill,
        ylims = ylims,
        ax_kws = None,
        savefig_kws = None
    )
    if fontsize is not None:
        xtick_labels = ax.get_xticklabels()
        for label in xtick_labels:
            label.set_fontsize(fontsize)
        ax.yaxis.get_label().set_fontsize(fontsize)
        ax.tick_params(axis = 'y', which = 'major', labelsize = fontsize)
    return(ax)



### BAF manipulation ###


def remove_nan_baf(mv, how = "independent"):
    """Remove genes whose baf is NaN for each group.
    
    how : {'independent', 'intersect'}
        'independent': filter genes independently in each group.
        'intersect': only keep intersected genes after 'independent' mode.
    """
    mv = mv.copy()
    mv['idx_inter'] = range(mv.shape[0])
    mv['idx_intra'] = -1   # gene name can also be used as index if available.
    for i, name in enumerate(mv['X_name'].unique()):
        df = mv.loc[mv['X_name'] == name].copy()
        mv.loc[mv['X_name'] == name, 'idx_intra'] = range(df.shape[0])
    
    res = mv.loc[~np.isnan(mv['baf'])].copy()

    if how == 'intersect':
        mv = res.copy()
        all_idx = []
        for i, name in enumerate(mv['X_name'].unique()):
            df = mv[mv['X_name'] == name].copy()
            if i == 0:
                all_idx = df['idx_intra']
            else:
                all_idx = np.intersect1d(all_idx, df['idx_intra'])
        res = mv.loc[mv['idx_intra'].isin(all_idx)].copy()
      
    r = res[[c for c in res.columns if c not in ('idx_inter', 'idx_intra')]]
    return(r, res['idx_inter'].to_numpy())



### Metrics ###


def get_baf_scalar(a, b):
    s = a + b
    if s == 0:
        return np.nan
    else:
        return b / float(s)
    

def get_baf(A, B):
    res = []
    for a, b in zip(A, B):
        res.append(get_baf_scalar(a, b))
    return np.array(res)



def __get_gw_metrics(X, metrics = None):        
    all_metrics = ["mean", "var", "cv", "zero_prop", "sum"]
    if metrics is None:
        metrics = all_metrics
    for m in metrics:
        if m not in all_metrics:
            error("invalid metric '%s'." % m)
            raise ValueError
            
    X = sparse2array(X)
    
    res = []
    for m in metrics:
        if m == "mean":
            res.append(get_gw_mean(X))
        elif m == "var":
            res.append(get_gw_var(X))
        elif m == "cv":
            res.append(get_gw_cv(X))
        elif m == "zero_prop":
            res.append(get_gw_zero_prop(X))
        elif m == "sum":
            res.append(get_gw_sum(X))
        else:
            error("invalid metric '%s'." % m)
            raise ValueError

    res = pd.DataFrame(data = {m:v for m, v in zip(metrics, res)})
    return(res)



def __get_gw_metrics_group(X_lst, id_lst, X_names = None, metrics = None):
    if X_names is None:
        X_names = ["X" + str(i) for i in range(len(X_lst))]
    if len(X_lst) != len(X_names):
        error("length of 'X_lst' and 'X_names' should be the same!")
        raise ValueError
    
    result = None
    for idx, (X, ids, name) in enumerate(zip(X_lst, id_lst, X_names)):
        res = __get_gw_metrics(X, metrics = metrics)
        res["X_name"] = name
        res['feature'] = ids
        if idx == 0:
            result = res
        else:
            result = pd.concat([result, res], ignore_index = True)

    return(result)



def get_gw_metrics_group(adata_lst, id_lst, X_names, metrics = None):
    assert len(adata_lst) == len(X_names)
    mv_A = __get_gw_metrics_group(
        X_lst = [adata.layers["A"] for adata in adata_lst],
        id_lst = id_lst,
        X_names = X_names,
        metrics = metrics
    )
    mv_B = __get_gw_metrics_group(
        X_lst = [adata.layers["B"] for adata in adata_lst],
        id_lst = id_lst,
        X_names = X_names,
        metrics = metrics
    )
    mv_AB = __get_gw_metrics_group(
        X_lst = [adata.layers["A"] + adata.layers["B"] for adata in adata_lst],
        id_lst = id_lst,
        X_names = X_names,
        metrics = metrics
    )
    mv_ABU = __get_gw_metrics_group(
        X_lst = [adata.layers["A"] + adata.layers["B"] + adata.layers["U"] \
                 for adata in adata_lst],
        id_lst = id_lst,
        X_names = X_names,
        metrics = metrics
    )

    mv_A.columns = [c + '_A' for c in mv_A.columns]
    mv_B.columns = [c + '_B' for c in mv_B.columns]
    mv_AB.columns = [c + '_AB' for c in mv_AB.columns]
    mv_ABU.columns = [c + '_ABU' for c in mv_ABU.columns]
    
    mv = pd.concat([mv_A, mv_B, mv_AB, mv_ABU], axis = 1)
    mv["X_name"] = mv["X_name_A"]
    mv["feature"] = mv["feature_A"]
    mv["DP"] = mv["sum_AB"]
    mv["baf"] = get_baf(mv["sum_A"], mv["sum_B"])
    
    return(mv)



def get_gw_cv(X):
    """Get gene-wise coefficient of variation (CV)."""
    return np.std(X, axis = 0) / (np.mean(X, axis = 0) + 1e-8)

def get_gw_mean(X):
    """Get gene-wise mean."""
    # res = np.mean(X, axis = 0)
    # res = np.asarray(res).reshape(-1)
    return np.mean(X, axis = 0)

def get_gw_sum(X):
    """Get gene-wise sum."""
    # res = np.sum(X, axis = 0)
    # res = np.asarray(res).reshape(-1)
    return np.sum(X, axis = 0)

def get_gw_var(X):
    """Get gene-wise variance."""
    return np.var(X, axis = 0)

def get_gw_zero_prop(X):
    """Get gene-wise zero proportion."""
    return np.mean(X <= 0.0, axis = 0)
