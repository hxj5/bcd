# udf.py


import itertools
import numpy as np
import os
import pandas as pd


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
