# overlap.py - subset the prediction and truth label files given overlapping cells.


import gc
import numpy as np
import os
from logging import info
from ..utils.base import assert_e



def run_overlap(
    tool_list,
    tool_fn_list,
    truth_fn,
    out_dir, 
    out_prefix, 
    verbose = True
):
    """Subset the label files given overlapping cells.
    
    Parameters
    ----------    
    tool_list : list of Tool
        A list of tool-specific :class:`~.tool.Tool` objects.
    tool_fn_list : list of str
        A list of tool-specific label files containing tumor predictions.
    truth_fn : str
        A TSV file storing ground truth of tumor labels.
    out_dir : str
        The output folder.
    out_prefix : str
        The prefix to output files.
    verbose : bool, default True
        Whether to show detailed logging information.
        
    Returns
    -------        
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    return overlap_isec_cells(
        tool_list = tool_list, 
        tool_fn_list = tool_fn_list, 
        truth_fn = truth_fn,
        out_dir = out_dir, 
        out_prefix = out_prefix, 
        verbose = verbose
    )



def overlap_isec_cells(
    tool_list, 
    tool_fn_list, 
    truth_fn,
    out_dir, 
    out_prefix, 
    verbose = True
):
    """Subset the data by intersected cells only.
    
    Parameters
    ----------
    tool_list
    tool_fn_list
    truth_fn
    out_dir
    out_prefix
    verbose
        See :func:`run_overlap` for details.
        
    Returns
    -------        
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    # check args.
    assert len(tool_list) > 0
    assert len(tool_fn_list) == len(tool_list)
    for fn in tool_fn_list:
        assert_e(fn)
    assert_e(truth_fn)
    os.makedirs(out_dir, exist_ok = True)
    

    # get overlap (intersected) cells among tools.
    if verbose:
        info("get overlap (intersected) cells from %d tool files ..." % \
            len(tool_fn_list))

    ovp_cells = None
    for i, (tool, fn) in enumerate(zip(tool_list, tool_fn_list)):
        df = pd.read_csv(fn, delimiter = '\t')
        if i == 0:
            ovp_cells = df["barcode"]
        else:
            ovp_cells = np.intersect1d(df["barcode"], ovp_cells)
    fn = os.path.join(out_dir, "%s.tools.intersect.cells.tsv" % out_prefix)
    np.savetxt(fn, ovp_cells, fmt = "%s", delimiter = '\n')

    if verbose:
        info("tools: %d overlap cells." % len(ovp_cells))


    # further overlap with truth cells.
    if verbose:
        info("further overlap with truth cells ...")

    truth = pd.read_csv(fn, delimiter = '\t')
    ovp_cells = np.intersect1d(truth.obs["barcode"], ovp_cells)
    fn = os.path.join(out_dir, 
            "%s.tools_and_truth.intersect.cells.tsv" % out_prefix)
    np.savetxt(fn, ovp_cells, fmt = "%s", delimiter = '\n')

    if verbose:
        info("truth: %d overlap cells." % len(ovp_cells))


    # subset data given overlap cells.
    if verbose:
        info("subset data given overlap cells ...")

    out_tool_fn_list = []
    for tool, tool_fn in zip(tool_list, tool_fn_list):
        tid = tool.tid
        df = pd.read_csv(fn, delimiter = '\t')
        df.index = df['barcode']
        df = df.loc[ovp_cells].copy()
        
        fn = os.path.join(out_dir, "%s.%s.tsv" % (out_prefix, tid))
        df.to_csv(fn, sep = '\t', index = False)
        out_tool_fn_list.append(fn)
        
    truth = pd.read_csv(truth_fn, delimiter = '\t')
    truth.index = truth['barcode']
    truth = truth.loc[ovp_cells].copy()
        
    out_truth_fn = os.path.join(out_dir, "%s.truth.tsv" % (out_prefix, ))
    truth.to_csv(fn, sep = '\t', index = False)    


    res = dict(
        # out_tool_fn_list : list of str
        #   A list of output tool-specific label files, in the same order as 
        #   `tool_list` and `tool_fn_list`.
        out_tool_fn_list = out_tool_fn_list,
        
        # out_truth_fn : str
        #   Output subset truth label file.
        out_truth_fn = out_truth_fn
    )
    return(res)
