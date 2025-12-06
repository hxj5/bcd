# core.py


import os
import sys

from logging import info, error
from logging import warning as warn
from .steps.metric import run_metric
from .steps.overlap import run_overlap
from .steps.plot import run_plot
from .utils.base import assert_e
    
    
    
def bcd_core_pipeline(
    sid,
    tool_list,
    tool_fn_list,
    out_dir,
    truth_fn,
    overlap_how = 'isec',
    fig_dpi = 300,
    verbose = True
):
    """Pipeline.
    
    Parameters
    ----------
    sid : str
        Sample ID.
    tool_list : list of Tool
        A list of tool-specific :class:`~.tool.Tool` objects.
    tool_fn_list : list of str
        A list of tool-specific files storing predicted labels.
    out_dir : str
        The output folder.
    truth_fn : str
        An file storing ground truth labels.
    overlap_how
    fig_dpi
    verbose
        See :func:`~.main.bcd_main()` for details.
        
    Returns
    -------
    dict
        Results.
    """
    # check args.
    if len(tool_list) <= 0:
        info("no input tool data, skip all next steps ...")
        return(dict())

    assert len(tool_list) == len(tool_fn_list)
    for fn in tool_fn_list:
        assert_e(fn)
    os.makedirs(out_dir, exist_ok = True)

    step = 1

    
    # subset tool and truth labels given overlapping cells.
    info("subset tool and truth labels given overlapping cells ...")
    
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
        tool_fn_list = overlap_res["out_tool_fn_list"],
        truth_fn_list = overlap_res["out_truth_fn"],
        verbose = verbose
    )
    step += 1
    
    
    # plot metrics.
    info("plot metrics ...")
    
    res_dir = os.path.join(out_dir, "%d_plot" % step)
    os.makedirs(res_dir, exist_ok = True)
    plot_res = run_plot(
        sid,
        tool_list = tool_list,
        tool_fn_list = overlap_res["out_tool_fn_list"],
        truth_fn = overlap_res["out_truth_fn"],
        metric_fn = metric_res['out_fn'],
        out_dir = res_dir,
        out_prefix = sid,
        fig_dpi = fig_dpi,
        verbose = verbose
    )
    step += 1


    res = plot_res
    return(res)
