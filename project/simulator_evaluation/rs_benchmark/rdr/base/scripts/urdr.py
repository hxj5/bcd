# urdr.py


import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from ubase import (
    long2wide,        # udf
    get_log2FC,       # umath
    sparse2array,     # umatrix
    stat_LR1,
    stat_smooth       # usmooth
)



def stat_log2fc(
    data,
    metric = "mean", group = "group", gvars = None
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


def stat_stack_zeroprop(data, metric = "zero_prop", group = "group", 
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



### Linear Regression ###

def stat_LR(
    data,
    metric = "mean", group = "group", gvars = None, 
    log_scale = False
):
    """Statistics for linear regression.
    gvars (None or DataFrame).
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
        
    res = gvars[['x', 'y']].copy()
    res['a'] = np.nan
    res['b'] = np.nan
    res['r2'] = np.nan
    res['p_value'] = np.nan

    for i in range(res.shape[0]):
        cx = res['x'].loc[i]
        cy = res['y'].loc[i]
        a, b, r2, p_value = stat_LR1(
            x = df[cx].to_numpy(), 
            y = df[cy].to_numpy()
        )
        res.loc[i, 'a'] = a
        res.loc[i, 'b'] = b
        res.loc[i, 'r2'] = r2
        res.loc[i, 'p_value'] = p_value
    return(res)



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
    return np.std(X, axis = 0) / (np.mean(X, axis = 0) + 1e-10)

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
