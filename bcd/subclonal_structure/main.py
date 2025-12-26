# main.py


import os
import sys
import time

from logging import info, error
from logging import warning as warn
from .config import Config
from .core import bcd_core_pipeline
from .steps.predict import run_predict
from .steps.truth import run_truth
from .utils.base import assert_e
from .utils.xlog import init_logging
from ..app import APP, VERSION



def bcd_main(
    sid,
    tool_list,
    out_dir,
    truth_fn,
    n_cluster,
    overlap_how = 'isec',
    fig_dpi = 300,
    verbose = True
):
    """Main function.
    
    Parameters
    ----------
    sid : str
        Sample ID.
    tool_list : list of Tool
        A list of tool-specific :class:`~.tool.Tool` objects.
    out_dir : str
        The output folder.
    truth_fn : str
        A header-free file stroing the ground truth.
        Its first two columns should be:
        - `barcode` and `annotation`.
    n_cluster : int
        Number of clusters for tools that do not output clone labels.
    overlap_how : {"isec"}
        How to subset the tool matrices given the overlap cells.
        - "isec"
            Subset tool matrix by intersected cells only.
    fig_dpi : int, default 300
        Resolution of the plot.
    verbose : bool, default True
        Whether to show detailed logging information.
        
    Returns
    -------
    int
        The return code. 0 if success, negative otherwise.
    dict
        Results.
    """
    conf = Config()
    
    conf.sid = sid
    conf.tool_list = tool_list
    conf.out_dir = out_dir
    conf.truth_fn = truth_fn
    conf.n_cluster = n_cluster
    conf.overlap_how = overlap_how
    
    conf.fig_dpi = fig_dpi
    conf.verbose = verbose
    
    ret, res = bcd_run(conf)
    return((ret, res))

    
    
def bcd_run(conf):
    init_logging(stream = sys.stdout)

    ret = -1
    res = None

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)
    info("%s (VERSION %s)." % (APP, VERSION))

    try:
        res = bcd_core(conf)
    except ValueError as e:
        error(str(e))
        error("Running program failed.")
        error("Quiting ...")
        ret = -1
    else:
        info("All Done!")
        ret = 0
    finally:
        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        info("end time: %s" % time_str)
        info("time spent: %.2fs" % (end_time - start_time, ))

    return((ret, res))



def bcd_core(conf):
    bcd_init(conf)
    info("Configuration:")
    conf.show(fp = sys.stdout, prefix = "\t")
    
    
    pp_dir = os.path.join(conf.out_dir, "0_pp")
    os.makedirs(pp_dir, exist_ok = True)

    
    # predict subclonal structure.
    info("predict subclonal structure ...")
    
    #res_dir = os.path.join(pp_dir, "tools")
    res_dir = pp_dir
    os.makedirs(res_dir, exist_ok = True)
    predict_res = run_predict(
        tool_list = conf.tool_list, 
        out_dir = res_dir,
        k = conf.n_cluster,
        verbose = conf.verbose
    )
    
    
    # extract ground truth.
    info("extract ground truth ...")

    res_dir = os.path.join(pp_dir, "truth")
    os.makedirs(res_dir, exist_ok = True)
    truth_res = run_truth(
        truth_fn = conf.truth_fn, 
        out_fn = os.path.join(res_dir, "truth.tsv"),
        verbose = conf.verbose
    )
    
    
    # run core pipeline.
    info("run core pipeline ...")
        
    res = bcd_core_pipeline(
        sid = conf.sid,
        tool_list = conf.tool_list,
        tool_fn_list = predict_res["out_fns"],
        out_dir = conf.out_dir,
        truth_fn = truth_res['out_fn'],
        overlap_how = conf.overlap_how,
        fig_dpi = conf.fig_dpi,
        verbose = conf.verbose
    )
    return(res)



def bcd_init(conf):
    # check args.
    assert len(conf.tool_list) > 0

    os.makedirs(conf.out_dir, exist_ok = True)
    assert_e(conf.truth_fn)
