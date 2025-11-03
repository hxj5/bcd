# core.py


import os
import sys

from logging import info, error
from logging import warning as warn
from .steps.metric import run_metric
from .steps.overlap import run_overlap
from .steps.plot import run_plot
from .utils.base import assert_e
    
    
    
def bcd_cna_type(
    sid,
    cna_type,
    tool_list,
    tool_fn_list,
    out_dir,
    truth_fn,
    overlap_how = "isec-cells",
    max_n_cutoff = 1000,
    fig_width = 4.25,
    fig_height = 3.25,
    fig_dpi = 300,
    fig_dec = 3,
    fig_legend_xmin = 0.5,
    fig_legend_ymin = 0.12,
    verbose = True
):
    """Pipeline for one CNA type.
    
    Parameters
    ----------
    sid : str
        Sample ID.
    cna_type : str
        CNA types, one of {"gain", "loss", "loh"}.
    tool_list : list of Tool
        A list of tool-specific :class:`~.tool.Tool` objects.
    tool_fn_list : list of str
        A list of tool-specific adata files storing cell x gene matrices.
    out_dir : str
        The output folder.
    truth_fn : str
        An ".h5ad" File storing cell x gene ground truth binary matrix.
    overlap_how
    max_n_cutoff
    fig_width
    fig_height
    fig_dpi
    fig_dec
    fig_legend_xmin
    fig_legend_ymin
    verbose
        See :func:`~.main.bcd_main()` for details.
        
    Returns
    -------
    dict
        Results.
    """
    # check args.
    if len(tool_list) <= 0:
        info("%s: no input tool data, skip all next steps ..." % cna_type)
        return(dict())

    assert cna_type in ("gain", "loss", "loh")
    assert len(tool_list) == len(tool_fn_list)
    for fn in tool_fn_list:
        assert_e(fn)
    os.makedirs(out_dir, exist_ok = True)

    step = 1

    
    # subset the adata objects given overlapping cells and genes.
    info("subset the adata objects given overlapping cells and genes ...")
    
    res_dir = os.path.join(out_dir, "%d_overlap" % step)
    os.makedirs(res_dir, exist_ok = True)
    overlap_res = run_overlap(
        tool_list = tool_list, 
        tool_fn_list = tool_fn_list,
        truth_fn = truth_fn,
        overlap_how = overlap_how,
        out_dir = res_dir, 
        out_prefix = "overlap",
        verbose = verbose
    )
    step += 1
    
    
    # calculate metrics for each tool.
    info("calculate metrics for each tool ...")

    res_dir = os.path.join(out_dir, "%d_metric" % step)
    os.makedirs(res_dir, exist_ok = True)
    metric_res = run_metric(
        tool_list = tool_list,
        out_dir = res_dir,
        out_prefix = "metric",
        tool_fn_list = overlap_res["out_tool_fn_list"],
        truth_fn_list = overlap_res["out_truth_fn_list"],
        cna_type = cna_type,
        max_n_cutoff = max_n_cutoff,
        verbose = verbose
    )
    step += 1
    
    
    # plot metrics.
    info("plot metrics ...")
    
    res_dir = os.path.join(out_dir, "%d_plot" % step)
    os.makedirs(res_dir, exist_ok = True)
    plot_res = run_plot(
        sid = sid,
        cna_type = cna_type,
        out_dir = res_dir,
        out_prefix = sid,
        roc_fn = metric_res["roc_fn"],
        auroc_fn = metric_res["auroc_fn"],
        prc_fn = metric_res["prc_fn"],
        auprc_fn = metric_res["auprc_fn"],
        fig_width = fig_width,
        fig_height = fig_height,
        fig_dpi = fig_dpi,
        fig_dec = fig_dec,
        fig_legend_xmin = fig_legend_xmin,
        fig_legend_ymin = fig_legend_ymin,
        verbose = verbose
    )
    step += 1
 

    res = plot_res
    return(res)
