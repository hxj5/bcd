# uoutlier.py


import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scanpy as sc
import scipy as sp
import seaborn as sns
import statsmodels.api as sm
from scipy.stats import pearsonr, linregress
from scipy.stats import fisher_exact, false_discovery_control

from udf import long2wide
from umath import get_log2FC, sign, stat_LR
from uplot import savefig_wrapper



### Outlier-Related Functions ###

def detect_outliers(
    data, metric = "mean", group = "X_name", 
    log_scale = False, cn_ratio = 1.0, alpha = 0.05, min_log2FC = 1
):
    """Calculate the metrics of outliers in the pair-wise linear regressions.
    data (DataFrame).
    metric (str): column name in `data`. Column containing values for OLS.
    group (str): column name in `data`. That column contains group names.
    log_scale (bool): whether to transform `x` and `y` with log10(i+1).
    cn_ratio (float): copy number ratio of cancer cells to normal cells.
    alpha (float): FDR level for t-test.
    min_log2FC (float): minimum log2FC to be identified as outliers.
    """
    df = long2wide(data, columns = group, values = metric)
    df_lm = df
    if log_scale:
        df_lm = df.map(lambda x: np.log10(x + 1))
    groups = df.columns
    
    k = 0
    stat = None
    for i in range(df.shape[1] - 1):
        for j in range(i + 1, df.shape[1]):
            print("%d - %s vs. %s" % (k + 1, groups[j], groups[i]))
            t, d = stat_LR(x = df_lm.iloc[:, i], y = df_lm.iloc[:, j])
            x, y = df.iloc[:, i].copy(), df.iloc[:, j].copy()
            if groups[i] not in ("cs_cancer", "rs_cancer"):
                x *= cn_ratio
            if groups[j] not in ("cs_cancer", "rs_cancer"):
                y *= cn_ratio            
            fc = get_log2FC(x = x, y = y)
            res = pd.DataFrame(data = dict(
                feature_idx = range(t.shape[0]),
                x = groups[i],
                y = groups[j],
                value_x = df.iloc[:, i],
                value_y = df.iloc[:, j],
                adj_p = t[t.columns[-1]],
                cook_d = d[0],
                log2fc = fc
            ))
            if k == 0:
                stat = res
            else:
                stat = pd.concat([stat, res])
            k += 1
    stat.index = [str(i) for i in range(stat.shape[0])]
     
    idx = (stat["adj_p"] < alpha) & (np.abs(stat["log2fc"]) > min_log2FC)
    stat["is_outlier"] = 0
    stat.loc[idx, "is_outlier"] = 1
    return(stat, df)      # here df is wide-format.



def stat_heatmap_outliers(data, xvars, yvars):
    """Generate data for heatmap the outliers from the linear regressions.
    data (DataFrame): the long-format DataFrame containing outliers.
    xvars (list of str): x axis labels.
    yvars (list of str): y axis labels.
    """
    stat = data.groupby(['x', 'y'])['is_outlier'].agg(['sum', 'mean']).reset_index()
    stat.columns = ['x', 'y', 'count', 'frac']
    df = stat.pivot(index = 'x', columns = 'y', values = 'count')
    df = df.loc[xvars, yvars].transpose()
    df_frac = stat.pivot(index = 'x', columns = 'y', values = 'frac')
    df_frac = df_frac.loc[xvars, yvars].transpose()
    return(df, df_frac)
    

@savefig_wrapper
def heatmap_outliers(
    df, df_frac, 
    cmap, xrotation = 45, yrotation = 0,
    savefig_kws = None
):
    """Heatmap for the outliers from the linear regressions.
    df (DataFrame): wide-format DataFrame containing statistics for heatmap.
    df_frac (DataFrame).
    cmap (str): colormap name.
    """
    def __fmt_anno_n(i):
        if np.isnan(i): s = 'NA'
        else: s = '%d' % i
        return(s)
    
    def __fmt_anno_frac(i):
        if np.isnan(i): s = 'NA'
        elif i == 1: s = '100%'
        else: s = '%.2f%%' % (i*100, )
        s = '\n(%s)' % s
        return(s)
    
    anno = df.map(__fmt_anno_n) + df_frac.map(__fmt_anno_frac)
    mask = np.zeros_like(df)
    mask[np.triu_indices_from(mask, k = 1)] = True
    ax = sns.heatmap(df, cmap = cmap, mask = mask, square = True, 
                     annot = anno, fmt = 's')
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation = xrotation, rotation_mode = "anchor",
        ha = 'right'
    )
    ax.tick_params(axis = 'y', labelrotation = yrotation)
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    return(ax)



def __test_outlier_and_overlap(df):
    # table:
    #     outlier_overlap, outlier_non-overlap,
    #     non-outlier_overlap, non-outlier_non-overlap
    n11 = np.sum( (df["is_outlier"] == 1) & (df["is_overlap"] == 1) )
    n12 = np.sum( (df["is_outlier"] == 1) & (df["is_overlap"] == 0) )
    n21 = np.sum( (df["is_outlier"] == 0) & (df["is_overlap"] == 1) )
    n22 = np.sum( (df["is_outlier"] == 0) & (df["is_overlap"] == 0) )
    data = np.array([[n11, n12], [n21, n22]])
    odds_ratio, p_value = fisher_exact(data)
    return(pd.DataFrame(
        data = dict(
            n_outlier_overlap = n11, 
            n_outlier_nonoverlap = n12,
            n_nonoutlier_overlap = n21, 
            n_nonoutlier_nonoverlap = n22,            
            odds_ratio = odds_ratio, 
            p_value = p_value),
        index = [0]
    ))


def stat_outlier_and_overlap(data, cna_genes, ovp_genes):
    """Test a set of outliers for their associations with overlap genes.
    data (DataFrame).
    cna_genes (array-like of str): CNA or non-CNA region genes.
    ovp_genes (array-like of str): overlap genes.
    """
    df = data.copy()
    df["feature"] = cna_genes[df["feature_idx"]]
    df["is_overlap"] = df["feature"].map(lambda g: 1 if g in ovp_genes else 0)
    stat = df.groupby(['x', 'y']).apply(__test_outlier_and_overlap).reset_index()
    stat['adj_p'] = false_discovery_control(stat['p_value'])
    return(stat, df)


def stat_heatmap_outlier_and_overlap(data, xvars, yvars):
    """Generate data for heatmap the overlap-outliers genes.
    data (DataFrame).
    xvars (list of str): x axis labels.
    yvars (list of str): y axis labels.
    """
    df = data.copy()
    df["n_outlier"] = df["n_outlier_overlap"]  + df["n_outlier_nonoverlap"]
    df["frac"] = df["n_outlier_overlap"] / df["n_outlier"]
    df.loc[df["n_outlier"] == 0, "frac"] = np.nan
    
    df_n = df.pivot(index = 'x', columns = 'y', values = 'n_outlier_overlap')
    df_n = df_n.loc[xvars, yvars].transpose()
    
    df_frac = df.pivot(index = 'x', columns = 'y', values = 'frac')
    df_frac = df_frac.loc[xvars, yvars].transpose()
    
    df_p = df.pivot(index = 'x', columns = 'y', values = 'adj_p')
    df_p = df_p.loc[xvars, yvars].transpose()
    return(df_n, df_frac, df_p)


@savefig_wrapper
def heatmap_outlier_and_overlap(
    df_n, df_frac, df_p, 
    cmap, show_p = True, anno_font_size = 10, xrotation = 45, yrotation = 0,
    savefig_kws = None
):
    """Heatmap for the fraction of overlap genes in outliers.
    df_n (DataFrame): number of outlier-overlap genes.
    df_frac (DataFrame): fraction of outlier-overlap genes in all outliers.
    df_p (DataFrame): p-value of outlier-overlap association test.
    cmap (str): colormap name.
    show_p (bool): whether to show p-value.
    """
    def __fmt_anno_n(i):
        if np.isnan(i): s = 'NA'
        else: s = '%d' % i
        return(s)
    
    def __fmt_anno_frac(i):
        if np.isnan(i): s = 'NA'
        elif i == 1: s = '100%'
        else: s = '%.2f%%' % (i*100, )
        s = '\n(%s)' % s
        return(s)
    
    def __fmt_anno_p(i):
        if np.isnan(i): s = 'p=NA'
        elif i < 0.001: s = 'p<0.001'
        else: s = 'p=%.3f' % i
        s = '\n' + s
        return(s)
    
    anno = None
    if show_p:
        anno = df_n.map(__fmt_anno_n) + \
                df_frac.map(__fmt_anno_frac) + \
                df_p.map(__fmt_anno_p)
    else:
        anno = df_n.map(__fmt_anno_n) + \
                df_frac.map(__fmt_anno_frac)

    mask = np.zeros_like(df_n)
    mask[np.triu_indices_from(mask, k = 1)] = True
    ax = sns.heatmap(
        df_n, 
        cmap = cmap, mask = mask, square = True, 
        annot = anno, annot_kws = {"fontsize":anno_font_size}, fmt = 's')
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation = xrotation, rotation_mode = "anchor",
        ha = 'right'
    )
    ax.tick_params(axis = 'y', labelrotation = yrotation)
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    return(ax)



def __pairgrid_with_outliers(
    ax, 
    x, y, data, outliers, 
    pos_r, pos_fit, 
    color_fit, color_ref, color_outlier,
    point_size, line_width
):
    """Mapping function for the lower triangle part of PairGrid."""
    data = data.copy()
    data.index = range(data.shape[0])
    cx, cy = x, y
    x, y = data[cx], data[cy]
    
    flag = (outliers['x'] == cx) & (outliers['y'] == cy) & \
            (outliers['is_outlier'] == 1)
    idx = outliers.loc[flag, "feature_idx"].to_numpy()

    x1 = data[cx].iloc[np.setdiff1d(data.index, idx)]   # non-outlier points.
    y1 = data[cy].iloc[np.setdiff1d(data.index, idx)]
    x2 = data[cx].iloc[idx]     # outlier points.
    y2 = data[cy].iloc[idx]
    
    #ax.scatter(x1, y1, s = point_size, marker = '.')
    #ax.scatter(x2, y2, s = point_size, color = color_outlier, marker = '^')
    ax.plot(x1, y1, marker = '.', markersize = point_size, linestyle = 'none')
    ax.plot(x2, y2, marker = '^', markersize = point_size, linestyle = 'none', 
            color = color_outlier)
    
    x = sm.add_constant(x)
    model = sm.OLS(y, x).fit()
    a, b = model.params.iloc[1], model.params.iloc[0]
    r2 = model.rsquared
    p_value = model.f_pvalue
    
    x_vals = np.linspace(max(ax.get_xlim()[0], np.min(x)), 
                         min(ax.get_xlim()[1], np.max(x)), 100)
    ax.plot(x_vals, a*x_vals+b, '-', color = color_fit, linewidth = line_width)
    ax.plot(x_vals, x_vals, '--', color = color_ref, linewidth = line_width)
    ax.annotate(r"R$^2$=%.2f, P=%.2e" % (r2, p_value), xy = pos_r, 
                xycoords = ax.transAxes)
    ax.annotate("y=%.2fx%s%.2f" % (a, sign(b), abs(b)), xy = pos_fit, 
                xycoords = ax.transAxes, color = color_fit)
    
    
def __pairgrid_without_outliers(
    ax, 
    x, y, data, outliers, 
    pos_r, pos_fit, 
    color_fit, color_ref, color_outlier,
    point_size, line_width,
    is_tri_upper = False
):
    """Mapping function for the lower triangle part of PairGrid.
    is_tri_upper (bool): whether plot within the upper part of triangle.
    """
    data = data.copy()
    data.index = range(data.shape[0])
    cx, cy = x, y
    x, y = data[cx], data[cy]
    
    # note that the `outliers` only have meaningful values for lower triangle part.
    flag = None
    if is_tri_upper:
        flag = (outliers['x'] == cy) & (outliers['y'] == cx) & \
                (outliers['is_outlier'] == 1)
    else:
        flag = (outliers['x'] == cx) & (outliers['y'] == cy) & \
                (outliers['is_outlier'] == 1)
    idx = outliers.loc[flag, "feature_idx"].to_numpy()

    x1 = data[cx].iloc[np.setdiff1d(data.index, idx)]   # non-outlier points.
    y1 = data[cy].iloc[np.setdiff1d(data.index, idx)]

    #ax.scatter(x1, y1, s = point_size, marker = '.')
    ax.plot(x1, y1, marker = '.', markersize = point_size, linestyle = 'none')
    
    x, y = x1, y1
    x = sm.add_constant(x)
    model = sm.OLS(y, x).fit()
    a, b = model.params.iloc[1], model.params.iloc[0]
    r2 = model.rsquared
    p_value = model.f_pvalue
    
    x_vals = np.linspace(max(ax.get_xlim()[0], np.min(x)), 
                         min(ax.get_xlim()[1], np.max(x)), 100)
    ax.plot(x_vals, a*x_vals+b, '-', color = color_fit, linewidth = line_width)
    ax.plot(x_vals, x_vals, '--', color = color_ref, linewidth = line_width)
    ax.annotate(r"R$^2$=%.2f, P=%.2e" % (r2, p_value), xy = pos_r, 
                xycoords = ax.transAxes)
    ax.annotate("y=%.2fx%s%.2f" % (a, sign(b), abs(b)), xy = pos_fit, 
                xycoords = ax.transAxes, color = color_fit)



@savefig_wrapper
def mulax_with_outliers(
    data, outliers, pos_r, pos_fit, 
    metric = "mean", group = "X_name", gvars = None, 
    xlabel = None, ylabel = None,
    log_scale = False,
    color_fit = None, color_ref = None, color_outlier = None,
    point_size = 25, line_width = 1,
    title_font_size = 10,
    nrows = 1, ncols = 1,
    savefig_kws = None
):
    """Multiple axes for comparisons between a set of groups.
    gvars (None or DataFrame).
    See `pair_grid()` for details of params.
    """
    df = long2wide(data, columns = group, values = metric)
    if log_scale:
        df = df.map(lambda x: np.log10(x + 1))
        
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
        __pairgrid_with_outliers(
            ax = ax,
            x = gvars['x'].loc[i], y = gvars['y'].loc[i], data = df, 
            outliers = outliers,
            pos_r = pos_r, pos_fit = pos_fit, 
            color_fit = color_fit, color_ref = color_ref, 
            color_outlier = color_outlier,
            point_size = point_size, line_width = line_width
        )
        ax.set_title(gvars['y'].loc[i] + ' ~ ' + gvars['x'].loc[i], 
                     fontsize = title_font_size)
        k += 1
        
    g.set_xlabels(label = xlabel, clear_inner = True)
    g.set_ylabels(label = ylabel, clear_inner = True)
    return(g)



@savefig_wrapper
def mulax_without_outliers(
    data, outliers, pos_r, pos_fit, 
    metric = "mean", group = "X_name", gvars = None, 
    xlabel = None, ylabel = None,    
    log_scale = False,
    color_fit = None, color_ref = None, color_outlier = None,
    point_size = 25, line_width = 1,
    title_font_size = 10,
    nrows = 1, ncols = 1,
    savefig_kws = None
):
    """Multiple axes for comparisons between a set of groups.
    See `pair_grid()` for details of params.
    """
    df = long2wide(data, columns = group, values = metric)
    if log_scale:
        df = df.map(lambda x: np.log10(x + 1))

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
        __pairgrid_without_outliers(
            ax = ax,
            x = gvars['x'].loc[i], y = gvars['y'].loc[i], data = df, 
            outliers = outliers,
            pos_r = pos_r, pos_fit = pos_fit, 
            color_fit = color_fit, color_ref = color_ref, 
            color_outlier = color_outlier,
            point_size = point_size, line_width = line_width
        )
        ax.set_title(gvars['y'].loc[i] + ' ~ ' + gvars['x'].loc[i], 
                     fontsize = title_font_size)
        k += 1
        
    g.set_xlabels(label = xlabel, clear_inner = True)
    g.set_ylabels(label = ylabel, clear_inner = True)
    return(g)
