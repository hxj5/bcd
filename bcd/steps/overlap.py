# overlap.py - subset the adata objects given overlapping cells and genes.


import gc
import numpy as np
import os
from logging import info
from ..utils.base import assert_e
from ..utils.io import load_h5ad, save_h5ad



def run_overlap(
    args_list, 
    tool_fn_list, 
    truth_fn,
    overlap_how,
    out_dir, 
    out_prefix, 
    verbose = True
):
    """Subset the adata objects given overlapping cells and genes.
    
    Parameters
    ----------    
    args_list : list of ToolArgs
        A list of tool-specific :class:`~.args.ToolArgs` objects.
    tool_fn_list : list of str
        A list of tool-specific adata files storing cell x gene matrices.
    truth_fn : str
        An ".h5ad" File storing cell x gene ground truth binary matrix.
    overlap_how : {"isec-cells", "isec-both"}
        How to subset the tool matrices given the overlap cells and genes.
        "isec-cells"
            Subset tool matrix by intersected cells only.
        "isec-both"
            Subset tool matrix by intersected cells and genes.
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
    # check args.
    assert len(args_list) > 0
    assert len(tool_fn_list) == len(args_list)
    for fn in tool_fn_list:
        assert_e(fn)
    assert_e(truth_fn)
    assert overlap_how in ("isec-cells", "isec-both")
    os.makedirs(out_dir, exist_ok = True)
    

    # get overlap (intersected) cells and genes among tools.
    out_tool_fn_list = out_truth_fn_list = None
    
    if overlap_how == "isec-cells":
        res = overlap_isec_cells(
            args_list = args_list, 
            tool_fn_list = tool_fn_list, 
            truth_fn = truth_fn,
            out_dir = out_dir, 
            out_prefix = out_prefix, 
            verbose = verbose
        )
        out_tool_fn_list = res["out_tool_fn_list"]
        out_truth_fn_list = res["out_truth_fn_list"]
        
    elif overlap_how == "isec-both":
        res = overlap_isec_both(
            args_list = args_list, 
            tool_fn_list = tool_fn_list, 
            truth_fn = truth_fn,
            out_dir = out_dir, 
            out_prefix = out_prefix, 
            verbose = verbose
        )
        out_tool_fn_list = res["out_tool_fn_list"]
        out_truth_fn_list = [res["out_truth_fn"]] * len(out_tool_fn_list)

    else:
        raise ValueError()


    res = dict(
        # out_tool_fn_list : list of str
        #   A list of output tool-specific adata files, in the same order as 
        #   `args_list` and `tool_fn_list`.
        out_tool_fn_list = out_tool_fn_list,
        
        # out_truth_fn : list of str
        #   Output subset truth adata file for each tool, in the same order as
        #   `args_list` and `tool_fn_list`.
        out_truth_fn_list = out_truth_fn_list
    )
    return(res)



def overlap_isec_cells(
    args_list, 
    tool_fn_list, 
    truth_fn,
    out_dir, 
    out_prefix, 
    verbose = True
):
    """Subset the adata objects by intersected cells only.
    
    Parameters
    ----------
    args_list
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
    # get overlap (intersected) cells among tools.
    if verbose:
        info("get overlap (intersected) cells from %d tool files ..." % \
            len(tool_fn_list))

    ovp_cells = None
    for i, (args, fn) in enumerate(zip(args_list, tool_fn_list)):
        adata = load_h5ad(fn)
        if i == 0:
            ovp_cells = adata.obs["cell"]
        else:
            ovp_cells = np.intersect1d(adata.obs["cell"], ovp_cells)
        del adata
        gc.collect()

    if verbose:
        info("tools: %d overlap cells." % len(ovp_cells))

    fn = os.path.join(out_dir, "%s.tools.intersect.cells.tsv" % out_prefix)
    np.savetxt(fn, ovp_cells, fmt = "%s", delimiter = '\n')
    
    
    # further overlap with truth cells.
    if verbose:
        info("further overlap with truth cells ...")

    truth = load_h5ad(truth_fn)
    ovp_cells = np.intersect1d(truth.obs["cell"], ovp_cells)

    if verbose:
        info("truth: %d overlap cells." % len(ovp_cells))

    fn = os.path.join(out_dir, 
            "%s.tools_and_truth.intersect.cells.tsv" % out_prefix)
    np.savetxt(fn, ovp_cells, fmt = "%s", delimiter = '\n')

    
    # subset adata given overlap cells and genes.
    if verbose:
        info("subset adata given overlap cells and genes ...")

    out_tool_fn_list = []
    out_truth_fn_list = []
    for args, tool_fn in zip(args_list, tool_fn_list):
        tid = args.tid

        adata = load_h5ad(tool_fn)
        old_shape_tool = adata.shape
        adata.obs.index = adata.obs["cell"]
        adata.var.index = adata.var["gene"]
        
        truth = load_h5ad(truth_fn)
        old_shape_truth = truth.shape
        truth.obs.index = truth.obs["cell"]
        truth.var.index = truth.var["gene"]

        ovp_genes = np.intersect1d(adata.var["gene"], truth.var["gene"])
        fn = os.path.join(out_dir, 
                "%s.%s_and_truth.intersect.genes.tsv" % (out_prefix, tid))
        np.savetxt(fn, ovp_genes, fmt = "%s", delimiter = '\n') 
        
        adata = adata[ovp_cells, ovp_genes]
        adata = adata.copy()
        truth = truth[ovp_cells, ovp_genes]
        truth = truth.copy()
        
        if verbose:
            info("%s: shape from %s to %s; truth: from %s to %s;" % \
                 (tid, str(old_shape_tool), str(adata.shape), 
                  str(old_shape_truth), str(truth.shape)))
        
        fn = os.path.join(out_dir, "%s.%s.h5ad" % (out_prefix, tid))
        save_h5ad(adata, fn)
        out_tool_fn_list.append(fn)
        
        fn = os.path.join(out_dir, "%s.truth.for_%s.h5ad" % (out_prefix, tid))
        save_h5ad(truth, fn)
        out_truth_fn_list.append(fn)
        
        del adata
        del truth
        gc.collect()


    res = dict(
        # out_tool_fn_list : list of str
        #   A list of output tool-specific adata files, in the same order as 
        #   `args_list` and `tool_fn_list`.
        out_tool_fn_list = out_tool_fn_list,
        
        # out_truth_fn : list of str
        #   Output subset truth adata file for each tool, in the same order as
        #   `args_list` and `tool_fn_list`.
        out_truth_fn_list = out_truth_fn_list
    )
    return(res)
    
    
    
def overlap_isec_both(
    args_list, 
    tool_fn_list, 
    truth_fn,
    out_dir, 
    out_prefix, 
    verbose = True
):
    """Subset the adata objects by intersected cells and genes.
    
    Parameters
    ----------    
    args_list
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
    # get overlap (intersected) cells and genes among tools.
    if verbose:
        info("get overlap (intersected) cells and genes from %d tool files ..." % \
            len(tool_fn_list))

    ovp_cells = ovp_genes = None
    for i, (args, fn) in enumerate(zip(args_list, tool_fn_list)):
        adata = load_h5ad(fn)
        if i == 0:
            ovp_cells = adata.obs["cell"]
            ovp_genes = adata.var["gene"]
        else:
            ovp_cells = np.intersect1d(adata.obs["cell"], ovp_cells)
            ovp_genes = np.intersect1d(adata.var["gene"], ovp_genes)
        del adata
        gc.collect()

    if verbose:
        info("tools: %d overlap cells and %d overlop genes." % \
             (len(ovp_cells), len(ovp_genes)))

    fn = os.path.join(out_dir, "%s.tools.intersect.cells.tsv" % out_prefix)
    np.savetxt(fn, ovp_cells, fmt = "%s", delimiter = '\n')
    
    fn = os.path.join(out_dir, "%s.tools.intersect.genes.tsv" % out_prefix)
    np.savetxt(fn, ovp_genes, fmt = "%s", delimiter = '\n')
    
    
    # further overlap with truth cells and genes.
    if verbose:
        info("further overlap with truth cells and genes ...")

    truth = load_h5ad(truth_fn)
    ovp_cells = np.intersect1d(truth.obs["cell"], ovp_cells)
    ovp_genes = np.intersect1d(truth.var["gene"], ovp_genes)
    
    if verbose:
        info("truth: %d overlap cells and %d overlop genes." % \
             (len(ovp_cells), len(ovp_genes)))

    fn = os.path.join(out_dir, 
            "%s.tools_and_truth.intersect.cells.tsv" % out_prefix)
    np.savetxt(fn, ovp_cells, fmt = "%s", delimiter = '\n')
    
    fn = os.path.join(out_dir, 
            "%s.tools_and_truth.intersect.genes.tsv" % out_prefix)
    np.savetxt(fn, ovp_genes, fmt = "%s", delimiter = '\n')
    
    
    # subset adata given overlap cells and genes.
    if verbose:
        info("subset adata given overlap cells and genes ...")

    out_tool_fn_list = []
    for args, tool_fn in zip(args_list, tool_fn_list):
        tid = args.tid

        adata = load_h5ad(tool_fn)
        old_shape = adata.shape
        
        adata.obs.index = adata.obs["cell"]
        adata.var.index = adata.var["gene"]
            
        adata = adata[ovp_cells, ovp_genes]
        adata = adata.copy()
        
        if verbose:
            info("%s: shape from %s to %s." % \
                 (tid, str(old_shape), str(adata.shape)))
        
        out_fn = os.path.join(out_dir, "%s.%s.h5ad" % (out_prefix, tid))
        save_h5ad(adata, out_fn)
            
        out_tool_fn_list.append(out_fn)
        
        del adata
        gc.collect()


    old_shape = truth.shape
    truth.obs.index = truth.obs["cell"]
    truth.var.index = truth.var["gene"]
    truth = truth[ovp_cells, ovp_genes]
    truth = truth.copy()
    
    if verbose:
        info("truth: shape from %s to %s." % \
                 (str(old_shape), str(truth.shape)))
        
    out_truth_fn = os.path.join(out_dir, "%s.%s.h5ad" % (out_prefix, "truth"))
    save_h5ad(truth, out_truth_fn)


    res = dict(
        # out_tool_fn_list : list of str
        #   A list of output tool-specific adata files, in the same order as 
        #   `args_list` and `tool_fn_list`.
        out_tool_fn_list = out_tool_fn_list,
        
        # out_truth_fn : str
        #   Output subset truth adata file.
        out_truth_fn = out_truth_fn
    )
    return(res)
