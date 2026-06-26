# ubase.py


import itertools
import numpy as np
import os
import pandas as pd
import scipy.sparse
import statsmodels.api as sm
from logging import info
from logging import warning as warn
from multiprocessing import Pool
from sccnasim.xlib.xmatrix import sparse2array
from scipy.interpolate import splrep, splev
from scipy.ndimage import gaussian_filter1d
from scipy.sparse import issparse
from statsmodels.nonparametric.kernel_regression import KernelReg



###### uDataFrame ######

def expand_grid(x, y, columns = None):
    if columns is None:
        columns = ['x', 'y']
    assert len(columns) == 2
    df = pd.DataFrame(itertools.product(x, y), columns = columns)
    return(df)



def long2wide(data, columns, values, index = None, keep_order = True):
    """Convert the long-format DataFrame to wide-format.
    Note that other columns other than `columns`, `values`, and `index` will
    be discarded.
    data (DataFrame).
    columns (str): column name in `data`. Column containing group names.
    values (str): column name in `data`. Column containing values.
    index (None or str): column name in `data`. Column containing indexes.
        If `None`, indexes will be set based on `columns` values.
    keep_order (bool): whether to keep the order of wide-format columns as in
        that in `columns`.
    """
    groups_col = data[columns].unique()
    groups_row = None
    df = None
    if index is None:
        df = data[[columns, values]].copy()
        n = len(groups_col)
        m = df.shape[0]
        assert m % n == 0
        df["idx"] = np.tile(range(int(m / n)), n)
        df = df.pivot(index = "idx", columns = columns, values = values)
    else:
        groups_row = data[index].unique()
        df = data[[index, columns, values]].copy()
        df = df.pivot(index = index, columns = columns, values = values)
    if keep_order:
        if index is None:
            df = df[groups_col]
        else:
            df = df.loc[groups_row, groups_col]
    return(df)



###### uMath ######

def format_equation(a, b, sep = ""):
    """Return formatted string for equation y=ax+b."""
    s = None
    if a == 0:
        if b == 0:
            s = "y%s=%s0" % (sep, sep)
        else:
            s = "y%s=%s%.2f" % (sep, sep, b)
    else:
        if b == 0:
            s = "y%s=%s%.2fx" % (sep, sep, a)
        else:
            s = "y%s=%s%.2fx%s%s%s%.2f" % (sep, sep, a, sep, sign(b), sep, abs(b))
    return(s)



def format_pvalue(p, sep = ""):
    """Return formatted string for p value."""
    s = "P%s=%s%.3f" % (sep, sep, p)
    if p < 0.001:
        s = "P%s<%s0.001" % (sep, sep)
    return(s)



def get_log2FC(x, y):
    """Calculate log2 fold change of each observation.
    x, y (array-like).
    """
    fc = np.zeros_like(x)
    fc[np.logical_and(x == 0, y != 0)] = np.inf
    fc[np.logical_and(x != 0, y == 0)] = -np.inf
    idx = np.logical_and(x != 0, y != 0)
    if np.sum(idx) > 0:
        fc[idx] = np.log2(y[idx] / x[idx])
    return(fc)



def mtx_norm_by_libsize(X, s = None):
    X = sparse2array(X)
    if s is None:
        s = X.sum(axis = 1)
    assert X.shape[0] == s.shape[0]
    s = s.reshape((s.shape[0], 1))
    X = X / s
    return(X)



def sign(x):
    return "+" if x >= 0 else "-"



def stat_LR1(x, y):
    x = sm.add_constant(x)
    model = sm.OLS(y, x).fit()
    try:
        a, b = model.params.iloc[1], model.params.iloc[0]
    except:
        a, b = model.params[1], model.params[0]
    r2 = model.rsquared
    p_value = model.f_pvalue
    return(a, b, r2, p_value)


  
# ref:
# - https://www.numberanalytics.com/blog/using-cooks-distance-advanced-outlier-detection
# - https://stackoverflow.com/questions/10231206/can-scipy-stats-identify-and-mask-obvious-outliers
def stat_LR_outlier(x, y):
    """Calculate p-value of t-test and Cook's distance of each observation
    in linear regression.
    
    x, y (array-like).
    """
    x = sm.add_constant(x)
    model = sm.OLS(y, x).fit()
    t_tests = model.outlier_test(method = "fdr_bh")
    influence = model.get_influence()
    cooks_d = influence.cooks_distance
    return(t_tests, cooks_d)



###### uMatrix ######

def array2sparse(X, which):
    """Convert a numpy array to a sparse array or matrix.
    
    Latest scipy library (v1.15.2) recommend to use sparse array rather than
    sparse matrix.
    Therefore, this function first tries converting numpy array to sparse
    array, and then sparse matrix if previous trial failed.

    Parameters
    ----------
    X : numpy.ndarray
        The numpy array to be converted into sparse arrya or sparse matrix.
    which : {"coo", "csc", "csr"}
        Which type of sparse array or matrix to use?
        - "coo": A sparse array/matrix in COOrdinate format.
        - "csc": Compressed Sparse Column array/matrix.
        - "csr": Compressed Sparse Row array/matrix.
    
    Returns
    -------
    sparse array or sparse matrix.
        The convertted sparse array or sparse matrix.
    """
    assert which in ("coo", "csc", "csr")
    
    # anndata seems not accepting sparse array when saving the adata object
    # into file via write_h5ad ...
    #try:
    #    if which == "coo":
    #        X = scipy.sparse.coo_array(X)
    #    elif which == "csc":
    #        X = scipy.sparse.csc_array(X)
    #    elif which == "csr":
    #        X = scipy.sparse.csr_array(X)
    #except:
    #    if which == "coo":
    #        X = scipy.sparse.coo_matrix(X)
    #    elif which == "csc":
    #        X = scipy.sparse.csc_matrix(X)
    #    elif which == "csr":
    #        X = scipy.sparse.csr_matrix(X)
    if which == "coo":
        X = scipy.sparse.coo_matrix(X)
    elif which == "csc":
        X = scipy.sparse.csc_matrix(X)
    elif which == "csr":
        X = scipy.sparse.csr_matrix(X)    
    return(X)


    
def sparse2array(X):
    """Convert a sparse matrix to numpy array.

    Parameters
    ----------
    X
        A sparse matrix.
    
    Returns
    -------
    numpy.ndarray
        The converted matrix.
    """
    if issparse(X):
        try:
            X = X.A
        except:
            X = X.toarray()
    return(X)



def mtx2array1d(mtx):
    """Convert numpy matrix-like (internally 1d) into 1d array."""
    return(mtx.A1)



###### uSmooth ######

# Notes
# 1. the `geom_smooth()` in R ggplot2 package may give unreasonable result that
#    some y (zero_prop) values are negative after smoothing.
# 2. similarly, the `KernelReg` in statsmodels.nonparametric.kernel_regression
#    may also give some negative predicted y (zero_prop) values.


def stat_smooth(
    x, y, data,
    k,
    group = "group", 
    log_scale = (False, False),
    ncores = 1
):
    """Do smoothing on two metrics.
    x, y (str): column names in `data`. Metric names.
    data (DataFrame).
    group (str): column name in `data`.
    log_scale (tuple of bool): whether to transform `x` and `y` with 
        log10(i+1).
    """
    cx, cy = x, y
    data = data.copy()

    assert len(log_scale) == 2
    log_scale_x, log_scale_y = log_scale
    if log_scale_x:
        data[cx] = data[cx].map(lambda i: np.log10(i+1))
    if log_scale_y:
        data[cy] = data[cy].map(lambda i: np.log10(i+1))
        
    df = None
    groups = data[group].unique()
    for g in groups:
        #info("processing group '%s' ..." % g)
        print("processing group '%s' ..." % g)
        d = data[data[group] == g].copy()
        x, y = d[cx], d[cy]
        x_smooth, y_smooth = smooth_knn_gk(x, y, k = k, ncores = ncores)
        res = pd.DataFrame(data = {cx: x_smooth, cy: y_smooth})
        res[group] = g
        res = res[[group, cx, cy]]
        if df is None:
            df = res
        else:
            df = pd.concat([df, res], ignore_index = True)
    return(df)



def gaussian_kernel(x_i, x_j, bw):
    return gaussian_kernel_1d(x_i, x_j, bw)


def gaussian_kernel_1d(x_i, x_j, bw):
    return np.exp( -0.5 * ( (x_i-x_j)/bw )**2 )


def gaussian_kernel_2d(x_i, y_i, x_j, y_j, bw):
    return np.exp( -0.5 * ( (x_i-x_j)**2 + (y_i-y_j)**2 ) / bw**2 )



def __smooth_knn_gk1(x_i, y_i, df, k, bw, dist):
    if dist == 'xy':
        df['dist'] = (x_i - df['x'])**2 + (y_i - df['y'])**2
        df2 = df.sort_values(by = ['dist', 'x', 'y'])
        df2 = df2.iloc[:k].copy()
        bw2 = max( bw, np.sqrt( np.var(df2['x']) + np.var(df2['y']) ) )
        weights = gaussian_kernel_2d(x_i, y_i, df2['x'], df2['y'], bw2)
        normalized_weights = weights / np.sum(weights)
        y_new = np.sum(normalized_weights * df2['y'])
        return(y_new)
    elif dist == 'x':
        df['dist'] = np.abs(x_i - df['x'])
        df2 = df.sort_values(by = ['dist', 'x', 'y'])
        df2 = df2.iloc[:k].copy()
        bw2 = max( bw, np.std(df2['x']) )
        weights = gaussian_kernel_1d(x_i, df2['x'], bw2)
        normalized_weights = weights / np.sum(weights)
        y_new = np.sum(normalized_weights * df2['y'])
        return(y_new)
    else:
        raise ValueError


def smooth_knn_gk(x, y, k, x_pred = None, n_pred = None, bw = None, 
                  dist = "x", ncores = 1):
    """KNN smoothing with Gaussian kernel."""
    n = len(x)
    if k > n:
        warn("k (%d) is greater than n (%d)." % (k, n))
        k = n
    
    # sort x and y.
    df = pd.DataFrame(data = dict(x = x, y = y))
    df = df.sort_values(by = ['x', 'y'])
    x, y = df['x'].to_numpy(), df['y'].to_numpy()
    x_pred, y_pred = None, None
    
    pool = Pool(processes = ncores)
    if dist == "xy":
        if bw is None:
            bw = np.sqrt(np.var(x) + np.var(y))
        assert x_pred is None
        assert n_pred is None
        
        # scenario 1.
        x_pred = x
        results = []
        for x_i, y_i in zip(x, y):
            res = pool.apply_async(
                __smooth_knn_gk1,
                args = (x_i, y_i, df, k, bw, dist)
            )
            results.append(res)
        pool.close()
        pool.join()

        y_pred = [res.get() for res in results]

    elif dist == "x":
        if bw is None:
            bw = np.std(x)
        if x_pred is None:
            if n_pred is None:       # scenario 2.
                x_pred = x
            else:                    # scenario 3.
                x_pred = np.linspace(x.min(), x.max(), num = n_pred)
        else:                        # scenario 3.
            pass

        results = []
        for x_i in x_pred:
            res = pool.apply_async(
                __smooth_knn_gk1,
                args = (x_i, None, df, k, bw, dist)
            )
            results.append(res)
        pool.close()
        pool.join()

        y_pred = [res.get() for res in results]

    else:
        raise ValueError

    return x_pred, np.array(y_pred)
