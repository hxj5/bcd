# usmooth.py


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from logging import warning as warn
from scipy.interpolate import splrep, splev
from scipy.ndimage import gaussian_filter1d
from statsmodels.nonparametric.kernel_regression import KernelReg

from uplot import ax_format_wrapper, savefig_wrapper



# Notes
# 1. the `geom_smooth()` in R ggplot2 package may give unreasonable result that
#    some y (zero_prop) values are negative after smoothing.
# 2. similarly, the `KernelReg` in statsmodels.nonparametric.kernel_regression
#    may also give some negative predicted y (zero_prop) values.



def stat_smooth(
    x, y, data, 
    group = "X_name", 
    log_scale = (False, False),
    smooth_kws = None
):
    """Do smoothing on two metrics.
    x, y (str): column names in `data`. Metric names.
    data (DataFrame).
    group (str): column name in `data`.
    log_scale (tuple of bool): whether to transform `x` and `y` with 
        log10(i+1).
    smooth_kws (dict): options passed to smoothing function.
    """
    cx, cy = x, y
    data = data.copy()

    assert len(log_scale) == 2
    log_scale_x, log_scale_y = log_scale
    if log_scale_x:
        data[cx] = data[cx].map(lambda i: np.log10(i+1))
    if log_scale_y:
        data[cy] = data[cy].map(lambda i: np.log10(i+1))
        
    if smooth_kws is None:
        smooth_kws = dict()
    
    # when using smooth_knn_gk() for smoothing
    assert 'k' in smooth_kws
        
    df = None
    groups = data[group].unique()
    for g in groups:
        d = data[data[group] == g].copy()
        x, y = d[cx], d[cy]
        x_smooth, y_smooth = smooth_knn_gk(x, y, **smooth_kws)
        res = pd.DataFrame(data = {cx: x_smooth, cy: y_smooth})
        res[group] = g
        res = res[[group, cx, cy]]
        if df is None:
            df = res
        else:
            df = pd.concat([df, res], ignore_index = True)
    return(df)
            


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
    """Smooth line plot.
    x, y (str): column names in `data`. Metric names.
    data (DataFrame).
    group (str): column name in `data`.
    legend (bool): whether to show legend.
    legend_pos (tuple of float): passed to `bbox_to_anchor`.
    """
    ax = None
    if legend is False:
        ax = sns.lineplot(data, x = x, y = y, hue = group, legend = legend,
                          linewidth = line_width)
    else:
        ax = sns.lineplot(data, x = x, y = y, hue = group, 
                          linewidth = line_width)
        if legend_pos is not None:
            sns.move_legend(
                ax, "upper left", title = None, frameon = False,
                bbox_to_anchor = legend_pos
            )
    if ylims is not None:
        ax.set_ylim(*ylims)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return(ax)



def gaussian_kernel(x_i, x_j, bw):
    return gaussian_kernel_1d(x_i, x_j, bw)


def gaussian_kernel_1d(x_i, x_j, bw):
    return np.exp( -0.5 * ( (x_i-x_j)/bw )**2 )


def gaussian_kernel_2d(x_i, y_i, x_j, y_j, bw):
    return np.exp( -0.5 * ( (x_i-x_j)**2 + (y_i-y_j)**2 ) / bw**2 )


def smooth_gkr1(x, y, bw = None, x_pred = None, n_pred = None, dist = "x", 
                itself = False):
    """Gaussian kernel regression, implementation from scratch.
    
    There are three distinct scenarios:
    1. predict the new `y` value for each of the input points (`x` and `y`),
       distance calculated based on `x` and `y`.
    2. predict the new `y` value for each of the input points (`x` and `y`),
       distance calculated based on `x` only.
    3. predict the new `y` value for each of the new points (`x_pred` or 
       from `n_pred`), distance calculated based on `x` only.
      
    Parameters
    ----------
    x : array-like
        The input `x` values.
    y : array-like
        The input `y` values.
    bw : float or None
        The bandwidth of the Gaussian kernel.
        If None, it will be inferred from the data.
    x_pred : array-like or None
        The x values of the predicted points.
        If both `x_pred` and `n_pred` are None, use input points for 
        prediction.
    n_pred : int or None
        The number of points to predict.
        At least one of `x_pred` and `n_pred` should be None.
    dist : str
        How distance or kernel is calculated.
        - "xy": use x and y values to calculate kernel.
        - "x": use only x values to calculate kernel.
    itself : bool
        Whether to count the point itself when calculating distance.
    
    Returns
    -------
    tuple
        array-like
            The newly predicted x values.
        array-like
            The newly predicted y values.
    """
    x_pred, y_pred = None, None
    
    # sort x and y.
    df = pd.DataFrame(data = dict(x = x, y = y))
    df = df.sort_values(by = ['x', 'y'])
    x, y = df['x'].to_numpy(), df['y'].to_numpy()
    
    if dist == "xy":
        if bw is None:
            bw = np.sqrt(np.var(x) + np.var(y))
        assert x_pred is None
        assert n_pred is None
        
        # scenario 1.
        x_pred = x
        y_pred = []
        for x_i, y_i in zip(x, y):
            weights = gaussian_kernel_2d(x_i, y_i, x, y, bw)
            normalized_weights = weights / np.sum(weights)
            y_i = np.sum(normalized_weights * y)
            y_pred.append(y_i)

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
        
        y_pred = []
        for x_i in x_pred:
            weights = gaussian_kernel_1d(x_i, x, bw)
            normalized_weights = weights / np.sum(weights)
            y_i = np.sum(normalized_weights * y)
            y_pred.append(y_i)
    
    else:
        raise ValueError

    return x_pred, np.array(y_pred)



def smooth_gkr2(x, y, x_pred = None, n_pred = None, **kws):
    """Use KernelReg() from statsmodels."""
    if x_pred is None:
        if n_pred is None:
            x_pred = np.sort(x)
        else:
            x_pred = np.linspace(x.min(), x.max(), num = n_pred)
    if 'var_type' not in kws:
        kws['var_type'] = 'c'
    kr = KernelReg(y, x, **kws)
    y_pred, _ = kr.fit(x_pred)
    return x_pred, np.array(y_pred)



def smooth_gkr3(x, y, sigma, **kws):
    """Implement with gaussian_filter1d() on x and y axis independently."""
    # ref: https://stackoverflow.com/questions/32900854/how-to-smooth-a-line-using-gaussian-kde-kernel-in-python-setting-a-bandwidth
    df = pd.DataFrame(data = dict(x = x, y = y))
    df = df.sort_values(by = ['x', 'y'])
    x, y = df['x'], df['y']
    x_smoothed = gaussian_filter1d(x, sigma, **kws)
    y_smoothed = gaussian_filter1d(y, sigma, **kws)
    return x_smoothed, y_smoothed



def smooth_knn_avg(x, y, k, x_pred = None, n_pred = None, dist = "x"):
    """KNN smoothing with simple average."""
    n = len(x)
    if k > n:
        warn("k (%d) is greater than n (%d)." % (k, n))
        k = n

    x_pred, y_pred = None, None
    
    # sort x and y.
    df = pd.DataFrame(data = dict(x = x, y = y))
    df = df.sort_values(by = ['x', 'y'])
    x, y = df['x'].to_numpy(), df['y'].to_numpy()
    
    if dist == "xy":
        assert x_pred is None
        assert n_pred is None
        
        # scenario 1.
        x_pred = x
        y_pred = []
        for x_i, y_i in zip(x, y):
            df['dist'] = (x_i - df['x'])**2 + (y_i - df['y'])**2
            df2 = df.sort_values(by = ['dist', 'y'])
            df2 = df2.iloc[:k].copy()
            y_i = np.mean(df2['y'])
            y_pred.append(y_i)

    elif dist == "x":
        if x_pred is None:
            if n_pred is None:       # scenario 2.
                x_pred = x
            else:                    # scenario 3.
                x_pred = np.linspace(x.min(), x.max(), num = n_pred)
        else:                        # scenario 3.
            pass
        
        y_pred = []
        for x_i in x_pred:
            df['dist'] = np.abs(x_i - df['x'])
            df2 = df.sort_values(by = ['dist', 'y'])
            df2 = df2.iloc[:k].copy()
            y_i = np.mean(df2['y'])
            y_pred.append(y_i)
    
    else:
        raise ValueError

    return x_pred, np.array(y_pred)



def smooth_knn_gk(x, y, k, x_pred = None, n_pred = None, bw = None, 
                  dist = "x"):
    """KNN smoothing with Gaussian kernel."""
    n = len(x)
    if k > n:
        warn("k (%d) is greater than n (%d)." % (k, n))
        k = n

    x_pred, y_pred = None, None
    
    # sort x and y.
    df = pd.DataFrame(data = dict(x = x, y = y))
    df = df.sort_values(by = ['x', 'y'])
    x, y = df['x'].to_numpy(), df['y'].to_numpy()
    
    if dist == "xy":
        if bw is None:
            bw = np.sqrt(np.var(x) + np.var(y))
        assert x_pred is None
        assert n_pred is None
        
        # scenario 1.
        x_pred = x
        y_pred = []
        for x_i, y_i in zip(x, y):
            df['dist'] = (x_i - df['x'])**2 + (y_i - df['y'])**2
            df2 = df.sort_values(by = ['dist', 'x', 'y'])
            df2 = df2.iloc[:k].copy()
            bw2 = max( bw, np.sqrt( np.var(df2['x']) + np.var(df2['y']) ) )
            weights = gaussian_kernel_2d(x_i, y_i, df2['x'], df2['y'], bw2)
            normalized_weights = weights / np.sum(weights)
            y_i = np.sum(normalized_weights * df2['y'])
            y_pred.append(y_i)

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
        
        y_pred = []
        for x_i in x_pred:
            df['dist'] = np.abs(x_i - df['x'])
            df2 = df.sort_values(by = ['dist', 'x', 'y'])
            df2 = df2.iloc[:k].copy()
            bw2 = max( bw, np.std(df2['x']) )
            weights = gaussian_kernel_1d(x_i, df2['x'], bw2)
            normalized_weights = weights / np.sum(weights)
            y_i = np.sum(normalized_weights * df2['y'])
            y_pred.append(y_i)
    
    else:
        raise ValueError

    return x_pred, np.array(y_pred)



def smooth_lowess(x, y, **kws):
    smoothed = sm.nonparametric.lowess(exog = x, endog = y, **kws)
    x_smoothed, y_smoothed = smoothed[:, 0], smoothed[:, 1]
    return x_smoothed, y_smoothed



def smooth_spline(x, y, s = 0, k = 3, **kws):
    # ref: https://stackoverflow.com/questions/78054656/what-is-the-difference-between-the-various-spline-interpolators-from-scipy
    # Note that it could be extreamly slow when data is big.
    df = pd.DataFrame(data = dict(x = x, y = y))
    df = df.sort_values(by = ['x', 'y'])
    x, y = df['x'], df['y']
        
    # typically spline requires `x` is strictly increasing for 1-d data.
    # for `x` with dup values, use below functions.
    ts = np.linspace(0, 1, len(y))
    ts_new = np.linspace(0, 1, len(y) * 2)
    (t0, c0, k0) = splrep(ts, x, s = s, k = k, **kws)
    (t1, c1, k1) = splrep(ts, y, s = s, k = k, **kws)
    x_smoothed = splev(ts_new, (t0, c0, k0))
    y_smoothed = splev(ts_new, (t1, c1, k1))
    return x_smoothed, y_smoothed



def smooth_plot_debug(x, y, x_smoothed, y_smoothed, title = None):
    ax = plt.gca()
    ax.scatter(x, y, label = 'Original', alpha = 0.7)
    ax.plot(x_smoothed, y_smoothed, color = 'red', linewidth = 2, 
            label = 'Smoothed')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if title is not None:
        ax.set_title(title)
    ax.legend()
    ax.grid(True)
    return(ax)
