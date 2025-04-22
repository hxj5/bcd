# main.py


import os
import sys
import time

from logging import info, error
from logging import warning as warn
from .app import APP, VERSION
from .config import Config
from .core import bcd_cna_type
from .steps.extract import run_extract
from .steps.truth import run_truth
from .utils.base import assert_e
from .utils.xlog import init_logging



def bcd_main(
    sid,
    args_list,
    out_dir,
    truth_fn,
    cell_anno_fn,
    gene_anno_fn,
    cna_type_list = None,
    max_n_cutoff = 1000,
    fig_width = 4.25,
    fig_height = 3.25,
    fig_dpi = 300,
    fig_dec = 3,
    fig_legend_xmin = 0.5,
    fig_legend_ymin = 0.12,
    verbose = True
):
    """Main function.
    
    Parameters
    ----------
    sid : str
        Sample ID.
    args_list : list of ToolArgs
        A list of tool-specific :class:`~.args.ToolArgs` objects.
    out_dir : str
        The output folder.
    truth_fn : str
        A header-free file stroing the ground truth.
        Its first five columns should be:
        - "chrom";
        - "start": 1-based and inclusive;
        - "end": 1-based and inclusive;
        - "clone": clone ID;
        - "cna_type": CNA type, should be in {"gain", "loss", "loh"}.
    cell_anno_fn : str
        File storing cell annotations.
        It is a header-free file whose first two columns are:
        - "cell": cell barcode;
        - "clone": clone ID;
    gene_anno_fn : str
        File storing gene annotations.
        It is a header-free file whose first four columns are:
        - "chrom";
        - "start": 1-based and inclusive;
        - "end": 1-based and inclusive;
        - "gene": gene name.
    cna_type_list : list of str or None, default None
        A list of CNA types.
        None means using all available CNA types, including "gain",
        "loss", and "loh".
    max_n_cutoff : int or None, default 1000
        Maximum number of cutoff values for calculating metrics.
        If None, use all unique values in tool matrix.
    fig_width : float, default 4.25
        Width of the plot in inch.
    fig_height : float, default 3.25
        Height of the plot in inch.
    fig_height : int, default 300
        Resolution of the plot.
    fig_dec : {3, 4}
        Number of decimal places for AUC.
    fig_legend_xmin : float, default 0.5
        The xmin position of legend.
    fig_legend_ymin : float, default 0.12
        The xmin position of legend.
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
    conf.args_list = args_list
    conf.out_dir = out_dir
    conf.truth_fn = truth_fn
    conf.cell_anno_fn = cell_anno_fn
    conf.gene_anno_fn = gene_anno_fn
    conf.cna_type_list = cna_type_list
    conf.max_n_cutoff = max_n_cutoff
    
    conf.fig_width = fig_width
    conf.fig_height = fig_height
    conf.fig_dpi = fig_dpi
    conf.fig_dec = fig_dec
    conf.fig_legend_xmin = fig_legend_xmin
    conf.fig_legend_ymin = fig_legend_ymin
    
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
    
    cna_type_dirs = []
    for cna_type in conf.cna_type_list:
        d = os.path.join(conf.out_dir, cna_type)
        os.makedirs(d, exist_ok = True)
        cna_type_dirs.append(d)

    
    # extract CNA expression or probability matrix into adata object.
    info("extract CNA expression or probability matrix into adata object ...")
    
    res_dir = os.path.join(pp_dir, "tools")
    os.makedirs(res_dir, exist_ok = True)
    extract_res = run_extract(
        args_list = conf.args_list, 
        out_dir = res_dir,
        out_prefix = "extract",
        cna_type_list = conf.cna_type_list,
        gene_anno_fn = conf.gene_anno_fn,
        verbose = conf.verbose
    )
    
    
    # extract CNA-type-specific ground truth.
    info("extract CNA-type-specific ground truth ...")

    res_dir = os.path.join(pp_dir, "truth")
    os.makedirs(res_dir, exist_ok = True)
    truth_res = run_truth(
        truth_fn = conf.truth_fn, 
        out_dir = res_dir, 
        out_prefix = "truth", 
        cna_type_list = conf.cna_type_list,
        cell_anno_fn = conf.cell_anno_fn,
        gene_anno_fn = conf.gene_anno_fn,
        verbose = conf.verbose
    )
    
    
    # run pipeline for each CNA type.
    info("run pipeline for each CNA type ...")
    
    cna_res = dict()
    for i, cna_type in enumerate(conf.cna_type_list):
        info("processing '%s' ..." % cna_type)
        
        args_list, tool_fn_list = [], []
        for args, fn in zip(conf.args_list, extract_res["out_fns"][cna_type]):
            if args.has_cna_type(cna_type) is True:
                args_list.append(args)
                tool_fn_list.append(fn)
        
        res = bcd_cna_type(
            sid = conf.sid,
            cna_type = cna_type,
            args_list = args_list,
            tool_fn_list = tool_fn_list,
            out_dir = cna_type_dirs[i],
            truth_fn = truth_res["out_fn_list"][i],
            max_n_cutoff = conf.max_n_cutoff,
            fig_width = conf.fig_width,
            fig_height = conf.fig_height,
            fig_dpi = conf.fig_dpi,
            fig_dec = conf.fig_dec,
            fig_legend_xmin = conf.fig_legend_xmin,
            fig_legend_ymin = conf.fig_legend_ymin,
            verbose = conf.verbose
        )
        cna_res[cna_type] = res
    

    res = cna_res
    return(res)



def bcd_init(conf):
    # check args.
    assert len(conf.args_list) > 0
    for args in conf.args_list:
        assert_e(args.obj_fn)

    os.makedirs(conf.out_dir, exist_ok = True)
    assert_e(conf.truth_fn)
    assert_e(conf.gene_anno_fn)
    
    if conf.cna_type_list is None:
        conf.cna_type_list = ("gain", "loss", "loh")
    else:
        assert len(conf.cna_type_list) > 0
        for cna_type in conf.cna_type_list:
            assert cna_type in ("gain", "loss", "loh")
