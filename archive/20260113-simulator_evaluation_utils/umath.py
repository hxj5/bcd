# umath.py


import numpy as np
import os
import pandas as pd
import statsmodels.api as sm
from sccnasim.xlib.xmatrix import sparse2array



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



def stat_LR(x, y):
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
