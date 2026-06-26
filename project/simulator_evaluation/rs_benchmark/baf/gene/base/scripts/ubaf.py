# ubaf.py


import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from ubase import sparse2array



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
    for i, name in enumerate(mv['group'].unique()):
        df = mv.loc[mv['group'] == name].copy()
        mv.loc[mv['group'] == name, 'idx_intra'] = range(df.shape[0])
    
    res = mv.loc[~np.isnan(mv['baf'])].copy()

    if how == 'intersect':
        mv = res.copy()
        all_idx = []
        for i, name in enumerate(mv['group'].unique()):
            df = mv[mv['group'] == name].copy()
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
    return np.std(X, axis = 0) / (np.mean(X, axis = 0) + 1e-10)

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
