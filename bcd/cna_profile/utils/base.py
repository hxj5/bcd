# base.py - basic utils.


import itertools
import numpy as np
import os
import pandas as pd
import subprocess


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



def is_function(x):
    """Test whether `x` is a function."""
    return hasattr(x, "__call__")


def is_scalar_numeric(x):
    """Test whether `x` is a scalar numeric value."""
    # ref: https://stackoverflow.com/questions/16807011/python-how-to-identify-if-a-variable-is-an-array-or-a-scalar
    return not hasattr(x, "__len__")


def is_vector(x):
    """Test whether `x` is a vector."""
    # ref: https://stackoverflow.com/questions/16807011/python-how-to-identify-if-a-variable-is-an-array-or-a-scalar
    return isinstance(x, (list, tuple, np.ndarray))


def is_file_empty(fn):
    """Test whether file is empty."""
    assert os.path.exists(fn)
    return(os.path.getsize(fn) <= 0)
