# cna_profile.py - CNA profile detection.

# Inputs
# * CalicoST - clonal gene copy number, spot-wise clone label and tumor prop.
# * CopyKAT - cell x gene CNA expression matrix;
# * InferCNV - cell x gene CNA expression matrix;
# * Numbat - cell x segment CNA probability;
# * XClone - cell x gene CNA probability;


import anndata as ad
import functools
import gc
import itertools
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy as sp
import seaborn as sns
import subprocess
import sys
import time
from logging import info, error
from logging import warning as warn
from sccnasim.xlib.xlog import init_logging
from sccnasim.xlib.xrange import format_chrom

from .app import APP, VERSION
#APP = "bcd"
#VERSION = "0.4.0"



############################################
#------------------ main ------------------#
############################################

def cna_profile_main(
    sid,
    tool_list,
    out_dir,
    truth_fn,
    cell_anno_fn,
    gene_anno_fn,
    cna_type_list = None,
    overlap_how = "isec-cells",
    max_n_cutoff = 1000,
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
    overlap_how : {"isec-cells", isec-both"}
        How to subset the tool matrices given the overlap cells and genes.
        - "isec-cells"
            Subset tool matrix by intersected cells only.
        - "isec-both"
            Subset tool matrix by intersected cells and genes.
    max_n_cutoff : int or None, default 1000
        Maximum number of cutoff values for calculating metrics.
        If None, use all unique values in tool matrix.
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
    conf.cell_anno_fn = cell_anno_fn
    conf.gene_anno_fn = gene_anno_fn
    conf.cna_type_list = cna_type_list
    conf.overlap_how = overlap_how
    conf.max_n_cutoff = max_n_cutoff
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
        tool_list = conf.tool_list,
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

        tool_list, tool_fn_list = [], []
        for tool, fn in zip(conf.tool_list, extract_res["out_fns"][cna_type]):
            if tool.has_cna_type(cna_type) is True:
                tool_list.append(tool)
                tool_fn_list.append(fn)

        res = bcd_cna_type(
            sid = conf.sid,
            cna_type = cna_type,
            tool_list = tool_list,
            tool_fn_list = tool_fn_list,
            out_dir = cna_type_dirs[i],
            truth_fn = truth_res["out_fn_list"][i],
            overlap_how = conf.overlap_how,
            max_n_cutoff = conf.max_n_cutoff,
            verbose = conf.verbose
        )
        cna_res[cna_type] = res


    res = cna_res
    return(res)



def bcd_init(conf):
    # check args.
    assert len(conf.tool_list) > 0

    os.makedirs(conf.out_dir, exist_ok = True)
    assert_e(conf.truth_fn)
    assert_e(conf.gene_anno_fn)

    if conf.cna_type_list is None:
        conf.cna_type_list = ("gain", "loss", "loh")
    else:
        assert len(conf.cna_type_list) > 0
        for cna_type in conf.cna_type_list:
            assert cna_type in ("gain", "loss", "loh")



##############################################
#------------------ config ------------------#
##############################################

class Config:
    def __init__(self):
        self.sid = None
        self.tool_list = None
        self.out_dir = None
        self.truth_fn = None
        self.cell_anno_fn = None
        self.gene_anno_fn = None
        self.cna_type_list = None
        self.overlap_how = "isec-cells"
        self.max_n_cutoff = 1000
        self.verbose = True


    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout

        s =  "%s\n" % prefix
        s += "%ssid = %s\n" % (prefix, self.sid)
        s += "%slen(tool_list) = %d\n" % (prefix, len(self.tool_list))
        s += "%stid list = '%s'\n" % (prefix, ", ".join([tool.display_name() for tool in self.tool_list]))
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%struth_fn = %s\n" % (prefix, self.truth_fn)
        s += "%scell_anno_fn = %s\n" % (prefix, self.cell_anno_fn)
        s += "%sgene_anno_fn = %s\n" % (prefix, self.gene_anno_fn)
        s += "%scna_type_list = %s\n" % (prefix, str(self.cna_type_list))
        s += "%soverlap_how = %s\n" % (prefix, self.overlap_how)
        s += "%smax_n_cutoff = %s\n" % (prefix, str(self.max_n_cutoff))
        s += "%sverbose = %s\n" % (prefix, str(self.verbose))
        s += "%s\n" % prefix

        fp.write(s)



############################################
#------------------ core ------------------#
############################################

def bcd_cna_type(
    sid,
    cna_type,
    tool_list,
    tool_fn_list,
    out_dir,
    truth_fn,
    overlap_how = "isec-cells",
    max_n_cutoff = 1000,
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
        verbose = verbose
    )
    step += 1


    res = plot_res
    return(res)



#####################################################
#------------------ steps.extract ------------------#
#####################################################

def _tool_file_id(tool):
    """Get unique filesystem identifier for tool (for output filenames)."""
    if tool.run_id is not None and str(tool.run_id).strip() != "":
        return "%s_%s" % (tool.tid.lower(), str(tool.run_id))
    return tool.tid.lower()


def run_extract(
    tool_list, out_dir, out_prefix,
    cna_type_list,
    gene_anno_fn,
    verbose = True
):
    """Extract CNA expression or probability matrix into adata object.

    Parameters
    ----------
    tool_list : list of Tools
        A list of tool-specific :class:`~..tools.Tool` objects.
    out_dir : str
        The output folder.
    out_prefix : str
        The prefix to output files.
    cna_type_list : list of str
        A list of CNA types, each in {"gain", "loss", "loh"}.
    gene_anno_fn : str
        File storing gene annotations.
    verbose : bool, default True
        Whether to show detailed logging information.

    Returns
    -------
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    # check args.
    assert len(tool_list) > 0
    os.makedirs(out_dir, exist_ok = True)
    assert len(cna_type_list) > 0
    assert_e(gene_anno_fn)


    out_fns = {cna_type:[] for cna_type in cna_type_list}
    for tool in tool_list:
        tid = tool.tid.lower()
        file_id = _tool_file_id(tool)
        info("extract matrix for '%s' ..." % tool.display_name())

        res_dir = os.path.join(out_dir, tid)
        os.makedirs(res_dir, exist_ok = True)

        if tid == "calicost":
            out_fn_list = [os.path.join(out_dir, "%s.%s.%s.h5ad" % \
                (out_prefix, file_id, cna_type)) for cna_type in cna_type_list]
            tool.extract(
                out_fn_list = out_fn_list,
                cna_type_list = cna_type_list,
                verbose = verbose
            )
            for cna_type, fn in zip(cna_type_list, out_fn_list):
                out_fns[cna_type].append(fn)

        elif tid == "copykat":
            fn = os.path.join(out_dir, "%s.%s.h5ad" % (out_prefix, file_id))
            tool.extract(
                out_fn = fn,
                verbose = verbose
            )
            for cna_type in cna_type_list:
                if tool.has_cna_type(cna_type):
                    out_fns[cna_type].append(fn)
                else:
                    out_fns[cna_type].append(None)

        elif tid == "infercnv":
            fn = os.path.join(out_dir, "%s.%s.h5ad" % (out_prefix, file_id))
            tool.extract(
                out_fn = fn,
                tmp_dir = res_dir,
                verbose = verbose
            )
            for cna_type in cna_type_list:
                if tool.has_cna_type(cna_type):
                    out_fns[cna_type].append(fn)
                else:
                    out_fns[cna_type].append(None)

        elif tid == "numbat":
            out_fn_list = [os.path.join(out_dir, "%s.%s.%s.h5ad" % \
                    (out_prefix, file_id, cna_type)) for cna_type in cna_type_list]
            tool.extract(
                out_fn_list = out_fn_list,
                cna_type_list = cna_type_list,
                gene_anno_fn = gene_anno_fn,
                tmp_dir = res_dir,
                verbose = verbose
            )
            for cna_type, fn in zip(cna_type_list, out_fn_list):
                out_fns[cna_type].append(fn)

        elif tid == "xclone":
            out_fn_list = [os.path.join(out_dir, "%s.%s.%s.h5ad" % \
                    (out_prefix, file_id, cna_type)) for cna_type in cna_type_list]
            tool.extract(
                out_fn_list = out_fn_list,
                cna_type_list = cna_type_list,
                verbose = verbose
            )
            for cna_type, fn in zip(cna_type_list, out_fn_list):
                out_fns[cna_type].append(fn)

        elif tid == "xclone_rdr":
            out_fn_list = [os.path.join(out_dir, "%s.%s.%s.h5ad" % \
                    (out_prefix, file_id, cna_type)) for cna_type in cna_type_list]
            tool.extract(
                out_fn_list = out_fn_list,
                cna_type_list = cna_type_list,
                verbose = verbose
            )
            for cna_type, fn in zip(cna_type_list, out_fn_list):
                out_fns[cna_type].append(fn)

        else:
            raise ValueError(f"Error: unknown tool id '{tid}'.")

    res = dict(
        # out_fns : dict of {str : list}
        #   Output CNA adata files.
        #   Each key is a CNA type, each value is a list of output adata files
        #   in the same order with `tool_list`.
        #   If one tool does not support some CNA type, then the adata file is
        #   set to `None`.
        out_fns = out_fns
    )

    return(res)



####################################################
#------------------ steps.metric ------------------#
####################################################

def run_metric(
    tool_list,
    out_dir,
    out_prefix,
    tool_fn_list,
    truth_fn_list,
    cna_type,
    max_n_cutoff = 1000,
    verbose = True
):
    """Main function of calculating metrics.

    Parameters
    ----------
    tool_list : list of Tool
        A list of tool-specific :class:`~.tool.Tool` objects.
    out_dir : str
        The output folder.
    out_prefix : str
        The prefix to output files.
    tool_fn_list : list of str
        A list of tool-specific files storing cell x gene expression or
        probability matrix.
    truth_fn_list : list of str
        A list of tool-specific files storing the ground truth cell x gene
        expression or probability matrices, matching order of `tool_fn_list`.
    cna_type : str
        CNA types, one of {"gain", "loss", "loh"}.
    max_n_cutoff : int or None, default 1000
        Maximum number of cutoff values for calculating metrics.
        If None, use all unique values in tool matrix.
    verbose : bool, default True
        Whether to show detailed logging information.

    Returns
    -------
    dict
        Results.
    """
    # check args.
    assert len(tool_list) == len(tool_fn_list)
    for fn in tool_fn_list:
        assert_e(fn)
    os.makedirs(out_dir, exist_ok = True)
    assert len(truth_fn_list) == len(tool_fn_list)
    for fn in truth_fn_list:
        assert_e(fn)
    assert cna_type in ("gain", "loss", "loh")


    # calculate metrics for each tool.
    if verbose:
        info("%s: calculate metrics for each tool ..." % cna_type)

    df_metric = None
    auroc_list = []
    auprc_list = []
    for i, (tool, tool_fn, truth_fn) in enumerate(
        zip(tool_list, tool_fn_list, truth_fn_list)):

        display_name = tool.display_name()
        if verbose:
            info("process %s ..." % display_name)

        adata = load_h5ad(tool_fn)
        tool_mtx = adata.X

        # NOTE: `tid` was undefined (bug). Use tool.tid instead.
        if tool.tid.lower() in ["infercnv", "copykat"] and cna_type == "loss":
            tool_mtx = tool_mtx * -1

        truth = load_h5ad(truth_fn)
        truth_mtx = truth.X

        res = calc_metric(
            tid = display_name,
            tool_mtx = tool_mtx,
            truth_mtx = truth_mtx,
            cna_type = cna_type,
            max_n_cutoff = max_n_cutoff,
            verbose = verbose
        )

        df = res["df"]
        df["tool"] = display_name
        if i == 0:
            df_metric = df.copy()
        else:
            df_metric = pd.concat([df_metric, df], ignore_index = True)
        auroc_list.append(res["AUROC"])
        auprc_list.append(res["AUPRC"])

        del adata
        del df
        del res
        gc.collect()


    # save files.
    if verbose:
        info("save metric result files ...")

    roc_fn = os.path.join(out_dir, "%s.%s.roc.tsv" % \
                             (out_prefix, cna_type))
    df_roc = df_metric[["tool", "cutoff", "tpr", "fpr"]]
    df_roc.to_csv(roc_fn, sep = "\t", index = False)


    prc_fn = os.path.join(out_dir, "%s.%s.prc.tsv" % \
                             (out_prefix, cna_type))
    df_prc = df_metric[["tool", "cutoff", "precision"]].copy()
    df_prc["recall"] = df_metric["tpr"]
    df_prc.to_csv(prc_fn, sep = "\t", index = False)


    df_auroc = pd.DataFrame(
        data = dict(
            tool = [tool.display_name() for tool in tool_list],
            AUROC = auroc_list
        ))
    auroc_fn = os.path.join(out_dir, "%s.%s.auroc.tsv" % \
                             (out_prefix, cna_type))
    df_auroc.to_csv(auroc_fn, sep = "\t", index = False)
    if verbose:
        info("%s: AUROC:\n\n%s\n" % (cna_type, str(df_auroc)))


    df_auprc = pd.DataFrame(
        data = dict(
            tool = [tool.display_name() for tool in tool_list],
            AUPRC = auprc_list
        ))
    auprc_fn = os.path.join(out_dir, "%s.%s.auprc.tsv" % \
                             (out_prefix, cna_type))
    df_auprc.to_csv(auprc_fn, sep = "\t", index = False)
    if verbose:
        info("%s: AUPRC:\n\n%s\n" % (cna_type, str(df_auprc)))


    res = dict(
        # roc_fn : str
        #   File storing DataFrame that contains ROC at various cutoffs.
        #   It has four columns "tool", "cutoff", tpr", "fpr".
        roc_fn = roc_fn,

        # prc_fn : str
        #   File storing DataFrame that contains PRC at various cutoffs.
        #   It has four columns "tool", "cutoff", "precision", "recall".
        prc_fn = prc_fn,

        # auroc_fn : str
        #   File storing DataFrame that contains AUROC of each tool.
        #   It has two columns "tool" and "AUROC".
        auroc_fn = auroc_fn,

        # auprc_fn : str
        #   File storing DataFrame that contains AUPRC of each tool.
        #   It has two columns "tool" and "AUPRC".
        auprc_fn = auprc_fn
    )
    return(res)



def calc_metric(
    tid,
    tool_mtx,
    truth_mtx,
    cna_type,
    max_n_cutoff = 1000,
    verbose = True
):
    """Calculate ROC and PRC.

    Parameters
    ----------
    tid : str
        Tool ID.
    tool_mtx : numpy.matrix-like or scipy.sparse matrix
        The tool's cell x gene expression or probability matrix.
    truth_mtx : numpy.matrix-like or scipy.sparse matrix
        The ground truth cell x gene expression or probability matrix.
    cna_type : str
        CNA types, one of {"gain", "loss", "loh"}.
    max_n_cutoff : int or None, default 1000
        Maximum number of cutoff values for `metric`.
        If None, use all unique values in `tool_mtx`.
    verbose : bool, default True
        Whether to show detailed logging information.

    Returns
    -------
    dict
        Results.
    """
    # check args.
    assert np.all(np.array(tool_mtx.shape) == np.array(truth_mtx.shape))
    assert cna_type in ("gain", "loss", "loh")

    if verbose:
        info("%s-%s: tool matrix %s and truth matrix %s." % \
             (cna_type, tid, str(tool_mtx.shape), str(truth_mtx.shape)))

    tool_vec = None
    if sp.sparse.issparse(tool_mtx):
        tool_vec = tool_mtx.toarray().flatten()
    else:
        tool_vec = tool_mtx.flatten()
    cutoff = np.unique(tool_vec)

    eps = (np.max(cutoff) - np.min(cutoff)) / len(cutoff)
    if eps <= 0:
        eps = 1e-6
    max_cutoff = np.max(cutoff) + eps
    min_cutoff = np.min(cutoff) - eps

    if max_n_cutoff is not None and len(cutoff) > max_n_cutoff - 2:
        cutoff = np.random.choice(cutoff, size = max_n_cutoff - 2, replace = False)
    cutoff = np.sort(np.unique(np.concatenate([cutoff, [min_cutoff, max_cutoff]])))

    if verbose:
        info("%s-%s: %d cutoffs." % (cna_type, tid, len(cutoff)))


    res = binary_metric(
        scores = tool_vec,
        labels = truth_mtx.toarray().flatten(),
        cutoff = cutoff,
        cut_direction = ">=",
        add_cut_inf = True
    )

    if verbose:
        info("%s-%s: AUROC = %f; AUPRC = %f;" % \
             (cna_type, tid, res["AUROC"], res["AUPRC"]))
    return(res)



# ref: the `binaryROC()` function in R package "cardelino".
# https://github.com/single-cell-genetics/cardelino/blob/3d2d9779877742dae8f3bc8644b5a0bee9973bd4/R/assessment.R#L192C1-L262C2
def binary_metric(
    scores,
    labels,
    cutoff = None,
    cut_direction = '>=',
    add_cut_inf = True
):
    """ROC and PRC curve for binary label prediction

    Parameters
    ----------
    scores : array-like of numeric
        Prediction score for each sample.
    labels : array-like of {0, 1}
        True labels for each sample.
    cutoff : array-like of numeric or None
        A vector of cutoffs; if None use all unique scores.
    cut_direction : {'>=', '>', '<=', '<'}
        A string to compare with cutoff.
    add_cut_inf : bool, default True
        If True, manually add a cutoff of `np.inf`.

    Returns
    -------
    dict
        Results.
    """
    # tpr = TP / (TP + FN)
    # fpr = FP / (FP + TN)
    # precision = TP / (TP + FP)
    # recall = tpr

    if cutoff is None:
        cutoff = np.sort(np.unique(scores))
    cutoff = np.sort(np.unique(np.concatenate([cutoff, [0, 1]])))

    n = len(cutoff)
    tpr = np.zeros(n)
    fpr = np.zeros(n)
    precision = np.zeros(n)

    idx = None
    for i in range(n):
        if cut_direction == "<=":
            idx = scores <= cutoff[i]
        elif cut_direction == "<":
            idx = scores < cutoff[i]
        elif cut_direction == ">":
            idx = scores > cutoff[i]
        else:
            idx = scores >= cutoff[i]

        if np.sum(idx) == 0:    # i.e., TP = FP = 0
            fpr[i] = 0.0
            tpr[i] = 0.0
            precision[i] = 1.0
        else:
            fpr[i] = np.sum(labels[idx] == 0) / np.sum(labels == 0)
            tpr[i] = np.sum(labels[idx] == 1) / np.sum(labels == 1)
            precision[i] = np.mean(labels[idx] == 1)
    #precision[logical_and(tpr == 0, np.isnan(precision))] = 1.0


    # Note,
    # For ROC and PRC plot, it is important to set the leftmost and rightmost
    # of the x axis, which is `fpr` (in ROC plot) and `tpr` (in PRC plot),
    # respectively.
    # - `fpr` and `tpr` decreases when cutoff increases.
    # - The `fpr` and `tpr` can always reach (1.0, 1.0) by setting an extreamly
    #   small cutoff;
    # - However, it cannot guarantee they can always reach (0.0, 0.0) unless
    #   manually set.
    if add_cut_inf:
        cutoff = np.concatenate([cutoff, [np.inf]])
        fpr = np.concatenate([fpr, [0.0]])
        tpr = np.concatenate([tpr, [0.0]])
        precision = np.concatenate([precision, [1.0]])
        n += 1

    auroc = auprc = 0.0
    for i in range(n - 1):
        auroc += (fpr[i] - fpr[i+1]) * (tpr[i] + tpr[i+1]) * 0.5
        auprc += (tpr[i] - tpr[i+1]) * (precision[i] + precision[i+1]) * 0.5
    auroc = auroc / (fpr[0] - fpr[-1])
    auprc = auprc / (tpr[0] - tpr[-1])

    df = pd.DataFrame(
        data = dict(
            cutoff = cutoff,
            tpr = tpr,
            fpr = fpr,
            precision = precision
        )
    )

    res = dict(
        # df : pandas.DataFrame
        #   A DataFrame containing AUROC and AUPRC at various cutoffs.
        df = df,

        # AUROC : float
        #   AUROC.
        AUROC = auroc,

        # AUPRC : float
        #   AUPRC.
        AUPRC = auprc
    )
    return(res)



#####################################################
#------------------ steps.overlap ------------------#
#####################################################

def run_overlap(
    tool_list,
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
    tool_list : list of Tool
        A list of tool-specific :class:`~.tool.Tool` objects.
    tool_fn_list : list of str
        A list of tool-specific adata files storing cell x gene matrices.
    truth_fn : str
        An ".h5ad" File storing cell x gene ground truth binary matrix.
    overlap_how : {"isec-cells", "isec-both"}
        How to subset the tool matrices given the overlap cells and genes.
        - "isec-cells"
            Subset tool matrix by intersected cells only.
        - "isec-both"
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
    assert len(tool_list) > 0
    assert len(tool_fn_list) == len(tool_list)
    for fn in tool_fn_list:
        assert_e(fn)
    assert_e(truth_fn)
    assert overlap_how in ("isec-cells", "isec-both")
    os.makedirs(out_dir, exist_ok = True)


    # get overlap (intersected) cells and genes among tools.
    out_tool_fn_list = out_truth_fn_list = None

    if overlap_how == "isec-cells":
        res = overlap_isec_cells(
            tool_list = tool_list,
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
            tool_list = tool_list,
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
        #   `tool_list` and `tool_fn_list`.
        out_tool_fn_list = out_tool_fn_list,

        # out_truth_fn : list of str
        #   Output subset truth adata file for each tool, in the same order as
        #   `tool_list` and `tool_fn_list`.
        out_truth_fn_list = out_truth_fn_list
    )
    return(res)



def overlap_isec_cells(
    tool_list,
    tool_fn_list,
    truth_fn,
    out_dir,
    out_prefix,
    verbose = True
):
    """Subset the adata objects by intersected cells only.

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
    # get overlap (intersected) cells among tools.
    if verbose:
        info("get overlap (intersected) cells from %d tool files ..." % \
            len(tool_fn_list))

    ovp_cells = None
    for i, (tool, fn) in enumerate(zip(tool_list, tool_fn_list)):
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
    for tool, tool_fn in zip(tool_list, tool_fn_list):
        file_id = _tool_file_id(tool)

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
                "%s.%s_and_truth.intersect.genes.tsv" % (out_prefix, file_id))
        np.savetxt(fn, ovp_genes, fmt = "%s", delimiter = '\n')

        adata = adata[ovp_cells, ovp_genes]
        adata = adata.copy()
        truth = truth[ovp_cells, ovp_genes]
        truth = truth.copy()

        if verbose:
            info("%s: shape from %s to %s; truth: from %s to %s;" % \
                 (tool.display_name(), str(old_shape_tool), str(adata.shape),
                  str(old_shape_truth), str(truth.shape)))

        fn = os.path.join(out_dir, "%s.%s.h5ad" % (out_prefix, file_id))
        save_h5ad(adata, fn)
        out_tool_fn_list.append(fn)

        fn = os.path.join(out_dir, "%s.truth.for_%s.h5ad" % (out_prefix, file_id))
        save_h5ad(truth, fn)
        out_truth_fn_list.append(fn)

        del adata
        del truth
        gc.collect()


    res = dict(
        # out_tool_fn_list : list of str
        #   A list of output tool-specific adata files, in the same order as
        #   `tool_list` and `tool_fn_list`.
        out_tool_fn_list = out_tool_fn_list,

        # out_truth_fn : list of str
        #   Output subset truth adata file for each tool, in the same order as
        #   `tool_list` and `tool_fn_list`.
        out_truth_fn_list = out_truth_fn_list
    )
    return(res)



def overlap_isec_both(
    tool_list,
    tool_fn_list,
    truth_fn,
    out_dir,
    out_prefix,
    verbose = True
):
    """Subset the adata objects by intersected cells and genes.

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
    # get overlap (intersected) cells and genes among tools.
    if verbose:
        info("get overlap (intersected) cells and genes from %d tool files ..." % \
            len(tool_fn_list))

    ovp_cells = ovp_genes = None
    for i, (tool, fn) in enumerate(zip(tool_list, tool_fn_list)):
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
    for tool, tool_fn in zip(tool_list, tool_fn_list):
        file_id = _tool_file_id(tool)

        adata = load_h5ad(tool_fn)
        old_shape = adata.shape

        adata.obs.index = adata.obs["cell"]
        adata.var.index = adata.var["gene"]

        adata = adata[ovp_cells, ovp_genes]
        adata = adata.copy()

        if verbose:
            info("%s: shape from %s to %s." % \
                 (tool.display_name(), str(old_shape), str(adata.shape)))

        out_fn = os.path.join(out_dir, "%s.%s.h5ad" % (out_prefix, file_id))
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
        #   `tool_list` and `tool_fn_list`.
        out_tool_fn_list = out_tool_fn_list,

        # out_truth_fn : str
        #   Output subset truth adata file.
        out_truth_fn = out_truth_fn
    )
    return(res)



##################################################
#------------------ steps.plot ------------------#
##################################################

def run_plot(
    sid,
    cna_type,
    out_dir,
    out_prefix,
    roc_fn,
    auroc_fn,
    prc_fn,
    auprc_fn,
    verbose = True
):
    """Main function of plotting ROC and PRC.

    Parameters
    ----------
    sid : str
        Sample ID.
    cna_type : str
        CNA types, one of {"gain", "loss", "loh"}.
    out_dir : str
        The output folder.
    out_prefix : str
        The prefix to output files.
    roc_fn : str
        File storing DataFrame that contains ROC at various cutoffs.
        It has four columns "tool", "cutoff", tpr", "fpr".
    auroc_fn : str
        File storing DataFrame that contains AUROC of each tool.
        It has two columns "tool" and "AUROC".
    prc_fn : str
        File storing DataFrame that contains PRC at various cutoffs.
        It has four columns "tool", "cutoff", "precision", "recall".
    auprc_fn : str
        File storing DataFrame that contains AUPRC of each tool.
        It has two columns "tool" and "AUPRC".
    verbose
        See :func:`~.main.bcd_main()` for details.

    Returns
    -------
    dict
        Results.
    """
    # check args.
    assert cna_type in ("gain", "loss", "loh")
    os.makedirs(out_dir, exist_ok = True)
    assert_e(roc_fn)
    assert_e(auroc_fn)
    assert_e(prc_fn)
    assert_e(auprc_fn)

    cna_name_map = {
        "gain" : "Copy Gain",
        "loss" : "Copy Loss",
        "loh" : "LOH"
    }

    out_fns = dict()
    for metric, df_fn, auc_fn in zip(("roc", "prc"), (roc_fn, prc_fn), \
                                     (auroc_fn, auprc_fn)):
        if verbose:
            info("plot %s ..." % metric.upper())

        out_fn = os.path.join(out_dir, "%s.%s.%s.jpg" % (out_prefix, cna_type, metric))
        plot_metric(
            metric = metric,
            df_fn = df_fn,
            auc_fn = auc_fn,
            out_fn = out_fn,
            title = "%s %s Curve for %s" % (sid, metric.upper(), cna_name_map[cna_type])
        )
        out_fns[metric] = out_fn


    res = dict(
        # out_fns : dict of {str : str}
        #   Output figures for each metric.
        #   Kyes are metrics, values are output figure files.
        out_fns = out_fns
    )
    return(res)



def plot_metric(
    metric,
    df_fn,
    auc_fn,
    out_fn,
    title = None,
    fig_width = 4.25,
    fig_height = 3.25,
    fig_dpi = 300,
    fig_dec = 3,
    fig_legend_xmin = 0.5,
    fig_legend_ymin = 0.12
):
    """Plot ROC or PRC curve.

    Parameters
    ----------
    metric : {"roc", "prc"}
        Evaluation metric.
    df_fn : str
        File storing DataFrame that contains ROC or PRC at various cutoffs.
        If ROC, it has four columns "tool", "cutoff", tpr", "fpr";
        otherwise, it has four columns "tool", "cutoff", "precision", "recall".
    auc_fn : str
        File storing DataFrame that contains AUROC or AUPRC of each tool.
        If AUROC, it has two columns "tool" and "AUROC";
        otherwise, it has two columns "tool" and "AUPRC".
    out_fn : str
        Output figure file.
    title : str or None
        Title of the plot figure. None means no title.
    fig_width
    fig_height
    fig_dpi
    fig_dec
    fig_legend_xmin
    fig_legend_ymin
        See :func:`~.main.bcd_main()` for details.

    Returns
    -------
    Void.
    """
    # check args.
    assert metric in ("roc", "prc")
    assert_e(df_fn)
    assert_e(auc_fn)


    df = pd.read_csv(df_fn, sep = "\t")
    auc = pd.read_csv(auc_fn, sep = "\t")

    if metric == "roc":
        for c in ("tool", "cutoff", "tpr", "fpr"):
            assert c in df.columns
    else:
        for c in ("tool", "cutoff", "precision", "recall"):
            assert c in df.columns

    if metric == "roc":
        for c in ("tool", "AUROC"):
            assert c in auc.columns
    else:
        for c in ("tool", "AUPRC"):
            assert c in auc.columns

    assert len(df["tool"].unique()) == len(auc["tool"].unique())
    assert np.all(np.sort(df["tool"].unique()) == np.sort(auc["tool"].unique()))


    # prepare plot data.
    df = df.merge(auc, on = "tool", how = "inner")

    if fig_dec == 4:
        if metric == "roc":
            df["legend"] = df.apply(lambda x: "%s: AUC=%.4f" % \
                                    (x["tool"], x["AUROC"]), axis = 1)
        else:
            df["legend"] = df.apply(lambda x: "%s: AUC=%.4f" % \
                                    (x["tool"], x["AUPRC"]), axis = 1)
    else:
        if metric == "roc":
            df["legend"] = df.apply(lambda x: "%s: AUC=%.3f" % \
                                    (x["tool"], x["AUROC"]), axis = 1)
        else:
            df["legend"] = df.apply(lambda x: "%s: AUC=%.3f" % \
                                    (x["tool"], x["AUPRC"]), axis = 1)


    # plot.
    title_fontsize = None
    axis_fontsize = None
    legend_fontsize = 8
    linewidth = 0.6

    fig = plt.figure(figsize = (fig_width, fig_height))

    if metric == "roc":
        sns.lineplot(data = df, x = 'fpr', y = 'tpr', hue = 'legend',
                        linewidth = linewidth, estimator = None, sort = False)
        plt.xlabel("False Positive Rate (1 - Specificity)", fontsize = axis_fontsize)
        plt.ylabel("True Positive Rate (Sensitivity)", fontsize = axis_fontsize)
    else:
        sns.lineplot(data = df, x = 'recall', y = 'precision', hue = 'legend',
                        linewidth = linewidth, estimator = None, sort = False)
        plt.xlabel("Recall", fontsize = axis_fontsize)
        plt.ylabel("Precision", fontsize = axis_fontsize)

    if title is not None:
        plt.title(title)

    if metric == "roc":
        plt.legend(
            loc = 'lower left',
            bbox_to_anchor=(fig_legend_xmin, fig_legend_ymin),
            fontsize = legend_fontsize
        )
    else:
        plt.legend(
            loc = 'lower right',
            bbox_to_anchor=(fig_legend_xmin, fig_legend_ymin),
            fontsize = legend_fontsize
        )

    plt.grid(False)
    plt.tight_layout()

    plt.savefig(out_fn, dpi = fig_dpi)
    plt.close(fig)



###################################################
#------------------ steps.truth ------------------#
###################################################

def run_truth(
    truth_fn,
    out_dir,
    out_prefix,
    cna_type_list,
    cell_anno_fn,
    gene_anno_fn,
    verbose = True
):
    """Format input ground truth table into cell x gene binary matrix.

    Parameters
    ----------
    truth_fn : str
        A header-free file stroing the ground truth.
        Its first five columns should be:
        - "chrom";
        - "start": 1-based and inclusive;
        - "end": 1-based and inclusive;
        - "clone": clone ID;
        - "cna_type": CNA type, should be in {"gain", "loss", "loh"}.
    out_dir : str
        The output folder.
    out_prefix : str
        The prefix to output files.
    cna_type_list : list of str
        A list of CNA types, each in {"gain", "loss", "loh"}.
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
    verbose : bool, default True
        Whether to show detailed logging information.

    Returns
    -------
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    # check args.
    assert_e(truth_fn)
    os.makedirs(out_dir, exist_ok = True)
    assert len(cna_type_list) > 0
    assert_e(cell_anno_fn)
    assert_e(gene_anno_fn)


    # merge truth profiles.
    if verbose:
        info("merge truth profiles ...")

    truth_fn_new = os.path.join(out_dir, "%s.merged.profiles.tsv" % out_prefix)
    ret, n_old, n_new = merge_truth_profiles(
        in_fn = truth_fn,
        out_fn = truth_fn_new,
        max_gap = 1
    )
    assert ret == 0
    if verbose:
        info("%d new CNA truth records merged from %d old ones." % \
             (n_new, n_old))


    # convert region-scale CNA truth profiles into gene-scale.
    if verbose:
        info("convert region-scale CNA truth profiles into gene-scale ...")

    df = load_truth(truth_fn_new)
    df["region"] = df.apply(
        lambda x: f"{x['chrom']}:{x['start']}-{x['end']}",
        axis = 1
    )
    df["region"] = df["region"].astype("object")
    df["is_cna"] = 1

    cell_anno = load_cell_anno(cell_anno_fn)
    assert len(cell_anno["cell"]) == len(cell_anno["cell"].unique())
    cell_anno["cell_index"] = range(cell_anno.shape[0])

    gene_anno = load_gene_anno(gene_anno_fn)
    assert len(gene_anno["gene"]) == len(gene_anno["gene"].unique())
    gene_anno["gene_index"] = range(gene_anno.shape[0])


    out_fn_list = []
    for cna_type in cna_type_list:
        if verbose:
            info("process '%s' ..." % cna_type)

        res_dir = os.path.join(out_dir, cna_type)
        os.makedirs(res_dir, exist_ok = True)

        c_df = df[df["cna_type"] == cna_type].copy()

        if verbose:
            s = "\n"
            s += "\t#unique regions = %d;\n" % len(c_df["region"].unique())
            s += "\t#unique clones = %d;\n" % len(c_df["clone"].unique())
            s += "\t#CNA records = %d;\n" % c_df["is_cna"].sum()

            info("%s: clone x region profile:\n%s" % (cna_type, s))


        # region-scale to gene-scale
        res = reg2gene(
            df = c_df,
            anno = gene_anno,
            no_dup = True,
            verbose = verbose
        )

        fn = os.path.join(res_dir, "overlap.mapping.tsv")
        res["overlap"].to_csv(fn, sep = "\t", index = False)

        c_df = res["df"].merge(cell_anno, on = "clone", how = "inner")
        assert c_df["is_cna"].isna().values.sum() == 0
        assert (c_df["is_cna"] != 1).sum() == 0

        fn = os.path.join(res_dir, "%s.truth.cell_x_gene.tsv" % cna_type)
        c_df.to_csv(fn, sep = "\t", index = False)

        if verbose:
            s = "\n"
            s += "\t#unique genes = %d;\n" % len(c_df["gene"].unique())
            s += "\t#unique cells = %d;\n" % len(c_df["cell"].unique())
            s += "\t#total records = %d;\n" % c_df.shape[0]
            s += "\t#CNA records = %d;\n" % c_df["is_cna"].sum()

            info("%s: cell x gene profile:\n%s" % (cna_type, s))


        # ground truth matrix should include all raw-data cells and genes.
        mtx = np.zeros((cell_anno.shape[0], gene_anno.shape[0]),
                       dtype = np.int8)

        assert c_df["is_cna"].isna().values.sum() == 0
        assert (c_df["is_cna"] != 1).sum() == 0

        row_idx = c_df["cell_index"].to_numpy()

        c_df = c_df.merge(
            gene_anno[["gene", "gene_index"]], on = "gene", how = "left")
        col_idx = c_df["gene_index"].to_numpy()

        mtx[row_idx, col_idx] = 1

        if verbose:
            info("%s: %d CNA records in matrix (shape = %s)." % \
                 (cna_type, mtx.sum(), str(mtx.shape)))


        # save truth matrix into adata file.
        adata = ad.AnnData(
            X = sp.sparse.csr_matrix(mtx),
            obs = pd.DataFrame(data = dict(cell = cell_anno["cell"])),
            var = pd.DataFrame(data = dict(gene = gene_anno["gene"]))
        )
        fn = os.path.join(out_dir, "%s.%s.h5ad" % (out_prefix, cna_type))
        save_h5ad(adata, fn)

        out_fn_list.append(fn)


    res = dict(
        # out_fn_list : list of str
        #   A list of CNA-type-specific ground truth ".h5ad" file.
        out_fn_list = out_fn_list
    )
    return(res)



def merge_truth_profiles(in_fn, out_fn, max_gap = 1):
    """Merge adjacent regions with the same CNA truth profiles.

    Merge adjacent regions with the same CNA types in each CNA clone.

    Parameters
    ----------
    in_fn : str
        Path to input file.
    out_fn : str
        Path to output file.
    max_gap : int, default 1
        The maximum gap length that is allowed between two adjacent regions.
        `1` for strict adjacence.

    Returns
    -------
    int
        The return code. 0 if success, negative if error.
    int
        Number of records before merging.
    int
        Number of records after merging.
    """
    sep = "\t"
    n_old, n_new = -1, -1


    # load data
    df = load_truth(in_fn, sep = sep)
    n_old = df.shape[0]

    dat = {}
    for i in range(df.shape[0]):
        rec = df.iloc[i, ]
        chrom, clone, cna_type = rec["chrom"], rec["clone"], rec["cna_type"]
        if clone not in dat:
            dat[clone] = {}
        if chrom not in dat[clone]:
            dat[clone][chrom] = {}
        if cna_type not in dat[clone][chrom]:
            dat[clone][chrom][cna_type] = []
        dat[clone][chrom][cna_type].append((rec["start"], rec["end"]))


    # merge (clone-specific) adjacent CNAs.
    for clone, cl_dat in dat.items():
        for chrom, ch_dat in cl_dat.items():
            for cna_type in ch_dat.keys():
                iv_list = sorted(
                    ch_dat[cna_type],
                    key = functools.cmp_to_key(cmp_two_intervals)
                )
                s1, e1 = iv_list[0]
                new_list = []
                for s2, e2 in iv_list[1:]:
                    if s2 <= e1 + max_gap:    # overlap adjacent region
                        e1 = max(e1, e2)
                    else:                     # otherwise
                        new_list.append((s1, e1))
                        s1, e1 = s2, e2
                new_list.append((s1, e1))
                ch_dat[cna_type] = new_list


    # check whether there are (strictly) overlapping regions with
    # distinct profiles.
    for clone, cl_dat in dat.items():
        for chrom, ch_dat in cl_dat.items():
            iv_list = []
            for cna_type in ch_dat.keys():
                iv_list.extend(
                    [(s, e, cna_type) for s, e in ch_dat[cna_type]])
            iv_list = sorted(
                iv_list,
                key = functools.cmp_to_key(cmp_two_intervals)
            )
            s1, e1 = iv_list[0][:2]
            for iv in iv_list[1:]:
                s2, e2 = iv[:2]
                if s2 <= e1:    # overlap adjacent region
                    error("distinct CNA profiles '%s', (%d, %d) and (%d, %d)." %
                        (chrom, s1, e1, s2, e2))
                    return((-5, n_old, n_new))
            cl_dat[chrom] = iv_list


    # save profile
    n_new = 0
    fp = open(out_fn, "w")
    for clone in sorted(dat.keys()):
        cl_dat = dat[clone]
        for chrom in sorted(cl_dat.keys()):
            ch_dat = cl_dat[chrom]
            for s, e, cna_type in ch_dat:
                fp.write("\t".join([chrom, str(s), str(e), \
                    str(clone), str(cna_type)]) + "\n")
                n_new += 1
    fp.close()
    return((0, n_old, n_new))



def cmp_two_intervals(x1, x2):
    """Compare two intervals.

    Parameters
    ----------
    x1 : list of (int, int)
        The begin and end coordinates of the first interval.
    x2 : list of (int, int)
        The begin and end coordinates of the second interval.

    Returns
    -------
    int
        The return code.
        0 if equal; positive if `x1` is larger; negative if `x2` is larger.
    """
    s1, e1 = x1[:2]
    s2, e2 = x2[:2]
    if s1 == s2:
        if e1 == e2:
            return(0)
        else:
            return e1 - e2
    else:
        return s1 - s2



##################################################
#------------------ tools.base ------------------#
##################################################

class Tool:
    def __init__(
        self,
        tid = None,
        has_gain = True, has_loss = True, has_loh = True,
        run_id = None
    ):
        self.tid = tid
        self.run_id = run_id
        self.__has_cna_type = {
            'gain': has_gain,
            'loss': has_loss,
            'loh': has_loh
        }


    def display_name(self):
        """Return display name for plots/legends. Includes run_id when present."""
        if self.run_id is not None and str(self.run_id).strip() != "":
            return "%s_%s" % (self.tid, self.run_id)
        return self.tid


    def has_cna_type(self, cna_type):
        if cna_type in self.__has_cna_type:
            return(self.__has_cna_type[cna_type])
        else:
            raise ValueError("invalid CNA type '%s'" % cna_type)



######################################################
#------------------ tools.calicost ------------------#
######################################################

class CalicoST(Tool):
    def __init__(self, cnv_fn, clone_fn, run_id = None):
        """Initialize CalicoST."""
        super().__init__(
            tid = "CalicoST",
            has_gain = True,
            has_loss = True,
            has_loh = True,
            run_id = run_id
        )
        self.cnv_fn = cnv_fn
        self.clone_fn = clone_fn


    def extract(self, out_fn_list, cna_type_list, verbose = False):
        """Extract CalicoST data and convert it to cell x gene probability
        matrices.

        Probabilities are scaled by tumor proportion for each cell.

        Parameters
        ----------
        out_fn_list : list of str
            Output ".h5ad" files storing the cell x gene matrix, each per
            CNA type.
        cna_type_list : list of str
            A list of CNA types, each in {"gain", "loss", "loh"}.
        verbose : bool, default False
            Whether to show detailed logging information.

        Returns
        -------
        Void.
        """
        return calicost_extract_cna_prob(
            cnv_fn = self.cnv_fn,
            clone_fn = self.clone_fn,
            out_fn_list = out_fn_list,
            cna_type_list = cna_type_list,
            verbose = verbose
        )



def calicost_extract_cna_prob(
    cnv_fn,
    clone_fn,
    out_fn_list,
    cna_type_list,
    verbose = False
):
    if verbose:
        info("Checking arguments...")

    assert_e(cnv_fn)
    assert_e(clone_fn)
    assert len(out_fn_list) > 0
    assert len(cna_type_list) == len(out_fn_list)
    for cna_type in cna_type_list:
        assert cna_type in ("gain", "loss", "loh")


    # Load input data.
    if verbose:
        info("Loading CalicoST data ...")

    cnv_df = pd.read_csv(cnv_fn, sep = '\t')
    clone_df = pd.read_csv(clone_fn, sep = '\t')

    assert 'gene' in cnv_df.columns, \
        "cnv_genelevel.tsv must contain 'gene' column!"
    for col in ['BARCODES', 'clone_label']:
        assert col in clone_df.columns, \
            "clone_labels.tsv must contain '%s' column!" % col
    if 'tumor_proportion' not in clone_df.columns:
        if verbose:
            warn("clone_labels.tsv has no 'tumor_proportion' column; assuming 1 for all cells.")
        clone_df = clone_df.copy()
        clone_df['tumor_proportion'] = 1.0


    # Handle duplicate barcodes and genes
    if clone_df['BARCODES'].duplicated().any():
        n_dup = clone_df['BARCODES'].duplicated().sum()
        warn(f"Found {n_dup} duplicate barcodes in clone_labels.tsv. " \
                "Keeping first occurrence.")
        clone_df = clone_df.drop_duplicates(
            subset = 'BARCODES', keep = 'first')

    if cnv_df['gene'].duplicated().any():
        n_dup = cnv_df['gene'].duplicated().sum()
        warn(f"Found {n_dup} duplicate genes in cnv_genelevel.tsv. "  \
                "Keeping first occurrence.")
        cnv_df = cnv_df.drop_duplicates(
            subset = 'gene', keep = 'first')


    # Init output data.
    genes = cnv_df['gene'].tolist()
    cells = clone_df['BARCODES'].tolist()
    clone_labels = clone_df['clone_label'].values
    tumor_proportions = clone_df['tumor_proportion'].values
    unique_clones = np.unique(clone_labels)


    obs_df = pd.DataFrame(data = dict(cell = cells))
    var_df = pd.DataFrame(data = dict(gene = genes))

    # Ensure unique indices
    # obs_df = obs_df.set_index('cell', verify_integrity = True)
    # var_df = var_df.set_index('gene', verify_integrity = True)


    # Process each CNA type
    for cna_type, out_fn in zip(cna_type_list, out_fn_list):
        if verbose:
            info(f"Processing CNA type '{cna_type}'...")

        mtx = np.zeros((len(cells), len(genes)), dtype = np.float32)
        for clone in unique_clones:
            cell_mask = clone_df['clone_label'] == clone
            cell_indices = np.where(cell_mask)[0]

            col_a = f'clone{clone} A'
            col_b = f'clone{clone} B'
            assert col_a in cnv_df.columns, \
                'clone{clone}: missing A copy number in cnv_df!'
            assert col_b in cnv_df.columns, \
                'clone{clone}: missing B copy number in cnv_df!'

            total_cn = cnv_df[col_a] + cnv_df[col_b]

            # Assign base probabilities based on CNA type
            if cna_type == "loss":         # Loss: A+B < 2
                prob = (total_cn < 2)
            elif cna_type == "loh":        # LOH: A+B = 2 and (A=0 or B=0)
                prob = ((total_cn == 2) & \
                        ((cnv_df[col_a] == 0) | (cnv_df[col_b] == 0)))
            elif cna_type == "gain":       # Gain: A+B > 2
                prob = (total_cn > 2)
            else:
                raise ValueError(f"Error: unknown CNA type '{cna_type}'.")
            prob = prob.astype(np.float32)

            # Scale probabilities by tumor proportion for each cell
            # in the clone
            # CHECK ME: P(spot belongs to specific clone) == tumor_prop?
            for cell_idx in cell_indices:
                mtx[cell_idx, :] = prob * tumor_proportions[cell_idx]

        adata = ad.AnnData(
            X = mtx,
            obs = obs_df,
            var = var_df
        )

        save_h5ad(adata, out_fn)
        if verbose:
            info(f"Saved adata shape = '%s' for CNA type '%s'." %  \
                 (str(adata.shape), cna_type))

        del adata
        gc.collect()

    return(out_fn_list)



#####################################################
#------------------ tools.copykat ------------------#
#####################################################

class CopyKAT(Tool):
    def __init__(self, expr_mtx_fn, run_id = None):
        super().__init__(
            tid = "CopyKat",
            has_gain = True,
            has_loss = True,
            has_loh = False,
            run_id = run_id
        )
        self.expr_mtx_fn = expr_mtx_fn


    def extract(self, out_fn, verbose = False):
        """
        Extract CNV data from CopyKAT output and save as AnnData object.

        Parameters
        ----------
        out_fn : str
            Output .h5ad file to save the matrix.
        verbose : bool, default False
            Whether to print verbose output.

        Returns
        -------
        dict
            {"mtx": cell x gene matrix (numpy array), "overlap": None}
        """
        return copykat_extract_cna_expression(
            expr_mtx_fn = self.expr_mtx_fn,
            out_fn = out_fn,
            verbose = verbose
        )



def copykat_extract_cna_expression(
    expr_mtx_fn,
    out_fn,
    verbose = False
):
    if verbose:
        info(f"Reading CopyKAT file: {expr_mtx_fn}")
    assert_e(expr_mtx_fn)

    mtx = pd.read_csv(expr_mtx_fn, sep = "\t", header = 0, dtype = str)
    mtx.index = mtx["hgnc_symbol"]
    mtx = mtx.iloc[:, 7:]              # Remove first 7 columns.
    mtx = mtx.T                        # cell x gene

    if verbose:
        info(f"CopyKAT matrix shape: {mtx.shape}")

    adata = ad.AnnData(
        X = mtx.values.astype(float),
        obs = pd.DataFrame(data = dict(cell = mtx.index)),
        var = pd.DataFrame(data = dict(gene = mtx.columns))
    )
    save_h5ad(adata, out_fn)

    if verbose:
        info(f"Saved AnnData to {out_fn}")

    return(out_fn)



######################################################
#------------------ tools.infercnv ------------------#
######################################################

class InferCNV(Tool):
    def __init__(self, obj_fn, run_id = None):
        """
        obj_fn : str
            File storing the inferCNV object.
            Typically using the "run.final.infercnv_obj".
        run_id : str or None, default None
            Optional run identifier for multiple runs of the same tool.
        """
        super().__init__(
            tid = "inferCNV",
            has_gain = True,
            has_loss = True,
            has_loh = False,
            run_id = run_id
        )
        self.obj_fn = obj_fn


    def extract(self, out_fn, tmp_dir, verbose = False):
        return infercnv_extract_cna_expression(
            obj_fn = self.obj_fn,
            out_fn = out_fn,
            tmp_dir = tmp_dir,
            verbose = verbose
        )



def infercnv_extract_cna_expression(obj_fn, out_fn, tmp_dir, verbose = False):
    """Extract inferCNV expression matrix and convert it to python object.

    Parameters
    ----------
    out_fn : str
        Output ".h5ad" file storing the cell x gene matrix.
    tmp_dir : str
        The folder to store temporary data.
    verbose : bool, default False
        Whether to show detailed logging information.

    Returns
    -------
    Void.
    """
    # check args.
    if verbose:
        info("check args ...")

    assert_e(obj_fn)
    os.makedirs(tmp_dir, exist_ok = True)


    # generate R scripts.
    if verbose:
        info("generate R scripts ...")

    cell_fn = os.path.join(tmp_dir, "barcodes.tsv")
    gene_fn = os.path.join(tmp_dir, "genes.tsv")
    mtx_fn = os.path.join(tmp_dir, "matrix.mtx")

    s = ""
    s += '''# extract cell x gene expression matrix from infercnv output.\n'''
    s += '''\n'''
    s += '''obj <- readRDS("%s")\n''' % obj_fn
    s += '''mtx <- obj@expr.data\n'''
    s += '''mtx <- t(mtx)         # cell x gene matrix\n'''
    s += '''\n'''
    s += '''write(\n'''
    s += '''    rownames(mtx),\n'''
    s += '''    file = "%s"\n''' % cell_fn
    s += ''')\n'''
    s += '''\n'''
    s += '''write(\n'''
    s += '''    colnames(mtx),\n'''
    s += '''    file = "%s"\n''' % gene_fn
    s += ''')\n'''
    s += '''\n'''
    s += '''write.table(\n'''
    s += '''    mtx,\n'''
    s += '''    file = "%s",\n''' % mtx_fn
    s += '''    row.names = FALSE,\n'''
    s += '''    col.names = FALSE\n'''
    s += ''')\n'''
    s += '''\n'''

    script_fn = os.path.join(tmp_dir, "extract_infercnv.R")
    with open(script_fn, "w") as fp:
        fp.write(s)


    # run the R script to save expression matrix into file.
    if verbose:
        info("run the R script to save expression matrix into file ...")

    exe_cmdline("Rscript %s" % script_fn)


    # load matrix into anndata.
    if verbose:
        info("load matrix into anndata and save into .h5ad file ...")

    barcodes = pd.read_csv(cell_fn, header = None)
    barcodes.columns = ["cell"]

    genes = pd.read_csv(gene_fn, header = None)
    genes.columns = ["gene"]

    mtx = np.loadtxt(mtx_fn)

    adata = ad.AnnData(
        X = mtx,
        obs = barcodes,
        var = genes
    )
    save_h5ad(adata, out_fn)

    if verbose:
        info("saved adata shape = %s." % str(adata.shape))

    return(out_fn)



####################################################
#------------------ tools.numbat ------------------#
####################################################

class Numbat(Tool):
    def __init__(self, joint_post_fn, mtx_how = 'expand', run_id = None):
        """
        joint_post_fn : str
            File storing the Numbat final results.
            Typically using the "joint_post_2.tsv".
        mtx_how : {"expand", "raw"}
            How to process the extracted Numbat matrix before overlap step.
            - "expand":
                expand the Numbat matrix to transcriptomics scale and fill value 0;
            - "raw":
                use the raw Numbat matrix.
        run_id : str or None, default None
            Optional run identifier for multiple runs of the same tool.
        """
        super().__init__(
            tid = "Numbat",
            has_gain = True,
            has_loss = True,
            has_loh = True,
            run_id = run_id
        )
        self.joint_post_fn = joint_post_fn
        self.mtx_how = mtx_how


    def extract(
        self,
        out_fn_list,
        cna_type_list,
        gene_anno_fn,
        tmp_dir,
        verbose = False
    ):
        """Extract Numbat probability matrices and convert them to python objects.

        Parameters
        ----------
        out_fn_list : list of str
            Output ".h5ad" files storing the cell x gene matrix, each per cna type.
        cna_type_list : list of str
            A list of CNA types, each in {"gain", "loss", "loh"}.
        gene_anno_fn : str
            File storing gene annotations.
        tmp_dir : str
            The folder to store temporary data.
        verbose : bool, default False
            Whether to show detailed logging information.

        Returns
        -------
        Void.
        """
        return numbat_extract_cna_prob(
            joint_post_fn = self.joint_post_fn,
            out_fn_list = out_fn_list,
            cna_type_list = cna_type_list,
            gene_anno_fn = gene_anno_fn,
            tmp_dir = tmp_dir,
            mtx_how = self.mtx_how,
            verbose = verbose
        )



def numbat_extract_cna_prob(
    joint_post_fn,
    out_fn_list,
    cna_type_list,
    gene_anno_fn,
    tmp_dir,
    mtx_how = 'expand',
    verbose = False
):
    # check args.
    if verbose:
        info("check args ...")

    assert_e(joint_post_fn)

    assert len(out_fn_list) > 0

    assert len(cna_type_list) == len(out_fn_list)
    for cna_type in cna_type_list:
        assert cna_type in ("gain", "loss", "loh")

    assert_e(gene_anno_fn)

    os.makedirs(tmp_dir, exist_ok = True)


    # load Numbat result.
    if verbose:
        info("load Numbat result ...")

    df = pd.read_csv(joint_post_fn, sep = '\t')
    df = df[["cell", "CHROM", "seg_start", "seg_end",
            "p_amp", "p_del", "p_loh", "p_bamp", "p_bdel"]]
    df.columns = ["cell", "chrom", "start", "end",
            "p_amp", "p_del", "p_loh", "p_bamp", "p_bdel"]

    df["chrom"] = df["chrom"].astype(str)
    df["region"] = df.apply(
        lambda x: f"{x['chrom']}:{x['start']}-{x['end']}",
        axis = 1
    )

    # sometimes Numbat outputs duplicate records (i.e., cell+region) when
    # there are multiple seg_labels, e.g., 1a_amp and 1a_loh.
    df = df.drop_duplicates(["cell", "region"], ignore_index = True)


    # get overlapping genes of each region.
    if verbose:
        info("get overlapping genes of each region ...")

    anno = load_gene_anno(gene_anno_fn)
    res = reg2gene(df, anno, verbose = verbose)

    fn = os.path.join(tmp_dir, "df.gene_scale.tsv")
    res["df"].to_csv(fn, sep = "\t", index = False)

    fn = os.path.join(tmp_dir, "overlap.mapping.tsv")
    res["overlap"].to_csv(fn, sep = "\t", index = False)


    df = res["df"][["cell", "gene",
                    "p_amp", "p_del", "p_loh", "p_bamp", "p_bdel"]]
    df = df.copy()

    cells = df["cell"].unique()
    genes = anno["gene"]
    df_ts = None
    if mtx_how == "expand":
        df_ts = pd.DataFrame(
            data = np.zeros((len(cells), len(genes)),
                            dtype = df["p_amp"].dtype),
            index = cells,
            columns = genes
        )

    for cna_type, out_fn in zip(cna_type_list, out_fn_list):
        if verbose:
            info("process cna_type '%s' ..." % cna_type)

        # calculate Numbat prob given `cna_type`.
        if cna_type == "gain":
            df["prob"] = df["p_amp"] + df["p_bamp"]
        elif cna_type == "loss":
            df["prob"] = df["p_del"] + df["p_bdel"]
        elif cna_type == "loh":
            df["prob"] = df["p_loh"]
        else:
            raise ValueError(f"Error: unknown cnv type '{cna_type}'.")

        # save gene-scale matrix into file.
        mtx = df.pivot(index = 'cell', columns = 'gene', values = 'prob')
        assert mtx.isna().values.sum() == 0
        X = mtx.to_numpy()

        if verbose:
            info("gene-scale matrix shape = %s." % str(X.shape))

        if mtx_how == "expand":
            df_tmp = df_ts.copy()
            df_tmp.loc[mtx.index, mtx.columns] = mtx
            mtx = df_tmp
            X = sp.sparse.csr_matrix(mtx.to_numpy())

        adata = ad.AnnData(
            X = X,
            obs = pd.DataFrame(data = dict(cell = mtx.index)),
            var = pd.DataFrame(data = dict(gene = mtx.columns))
        )
        save_h5ad(adata, out_fn)

        if verbose:
            info("saved adata shape = %s." % str(adata.shape))

        del adata
        gc.collect()

    return(out_fn_list)



####################################################
#------------------ tools.xclone ------------------#
####################################################

class XClone(Tool):
    def __init__(self, combine_fn, run_id = None):
        super().__init__(
            tid = "XClone",
            has_gain = True,
            has_loss = True,
            has_loh = True,
            run_id = run_id
        )
        self.combine_fn = combine_fn


    def extract(
        self,
        out_fn_list,
        cna_type_list,
        verbose = False
    ):
        """Extract XClone prob matrix and convert it to python object.

        Parameters
        ----------
        out_fn_list : list of str
            Output ".h5ad" files storing the cell x gene matrix, each per cna type.
        cna_type_list : list of str
            A list of CNA types, each in {"gain", "loss", "loh"}.
        verbose : bool, default False
            Whether to show detailed logging information.

        Returns
        -------
        Void.
        """
        return xclone_extract_cna_prob(
            combine_fn = self.combine_fn,
            out_fn_list = out_fn_list,
            cna_type_list = cna_type_list,
            verbose = verbose
        )



def xclone_extract_cna_prob(
    combine_fn,
    out_fn_list,
    cna_type_list,
    verbose = False
):
    if verbose:
        info("check args ...")

    assert_e(combine_fn)
    assert len(out_fn_list) > 0
    assert len(cna_type_list) == len(out_fn_list)

    for cna_type in cna_type_list:
        assert cna_type in ("gain", "loss", "loh")

    if verbose:
        info("load XClone object ...")

    xclone_adata = ad.read_h5ad(combine_fn)

    if verbose:
        info("XClone adata shape = %s." % str(xclone_adata.shape))

    prob = xclone_adata.layers['prob1_merge_refined']
    # prob_merge = np.stack([copy_loss, loh, copy_neutral, copy_gain], axis = -1)

    # Prepare obs and var DataFrames
    cells = xclone_adata.obs.index.tolist()
    genes = xclone_adata.var['GeneName'].tolist()

    obs_df = pd.DataFrame(data = dict(cell = cells))
    var_df = pd.DataFrame(data = dict(gene = genes))

    # iterate over cna types.
    for cna_type, out_fn in zip(cna_type_list, out_fn_list):
        if verbose:
            info("process cna_type '%s' ..." % cna_type)

        if cna_type == "loss":
            mtx = prob[:, :, 0]  # copy_loss
        elif cna_type == "loh":
            mtx = prob[:, :, 1]  # loh
        elif cna_type == "gain":
            mtx = prob[:, :, 3]  # copy_gain
        else:
            raise ValueError(f"Error: unknown cnv type '{cna_type}'.")

        adata = ad.AnnData(
            X = mtx,
            obs = obs_df,
            var = var_df
        )
        save_h5ad(adata, out_fn)

        if verbose:
            info("saved adata shape = %s." % str(adata.shape))

        del adata
        gc.collect()

    return(out_fn_list)



##################################################
#------------------ utils.base ------------------#
##################################################

def assert_e(path):
    """Assert file or folder exists, mimicking shell "test -e"."""
    assert path is not None and os.path.exists(path)



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



def is_file_empty(fn):
    """Test whether file is empty."""
    assert os.path.exists(fn)
    return(os.path.getsize(fn) <= 0)



####################################################
#------------------ utils.gscale ------------------#
####################################################

def get_overlap_genes(df, anno):
    """Get overlapping genes of regions.

    Parameters
    ----------
    df : pandas.DataFrame
        Region-scale data.
        It should contain at least four columns:
        - "chrom";
        - "start": 1-based and inclusive;
        - "end": 1-based and inclusive;
        - "region": region ID.
    anno : pandas.DataFrame
        Gene annotations.
        It should contain at least four columns:
        - "chrom";
        - "start": 1-based and inclusive;
        - "end": 1-based and inclusive;
        - "gene": gene name.

    Returns
    -------
    dict
        Overlapping results.
    """
    df["chrom"] = df["chrom"].map(format_chrom)
    anno["chrom"] = anno["chrom"].map(format_chrom)

    df = df.drop_duplicates("region", ignore_index = True)
    anno = anno.drop_duplicates("gene", ignore_index = True)

    overlap = None
    if df.shape[0] <= 0:
        overlap = pd.DataFrame(columns = ["region", "gene"], dtype = "object")
    else:
        overlap = df.groupby("region").apply(
            lambda x: anno.loc[
                (anno["chrom"] == x["chrom"].values[0]) &
                (anno["start"] <= x["end"].values[0]) &
                (anno["end"] >= x["start"].values[0]),
                ["chrom", "gene"]       # select two columns to force returning DataFrame instead of Series when there is only one region.
            ]).reset_index()
        overlap = overlap[["region", "gene"]].sort_values(by = "region")

    stat = overlap.groupby("gene").size().reset_index(name = "n")
    dup = stat[stat["n"] > 1]
    uniq = stat[stat["n"] == 1].merge(overlap, on = "gene", how = "left")
    uniq = uniq[["region", "gene"]].sort_values(by = "region")

    res = dict(
        # overlap : pandas.DataFrame
        #   The overlapping results. It contains two columns:
        #   - "region": region ID.
        #   - "gene": name of genes overlapping the region.
        overlap = overlap,

        # n_region : int
        #   Number of unique regions.
        n_region = df.shape[0],

        # n_region_overlap : int
        #   Number of unique regions that have overlapping genes.
        n_region_overlap = len(overlap["region"].unique()),

        # n_gene : int
        #   Number of unique genes.
        n_gene = anno.shape[0],

        # n_gene_overlap : int
        #   Number of unique genes that have overlapping regions.
        n_gene_overlap = len(overlap["gene"].unique()),

        # n_gene_dup : int
        #   Number of genes overlapping more than 1 regions.
        n_gene_dup = dup.shape[0],

        # overlap_uniq : pandas.DataFrame
        #   The overlapping results.
        #   Similar to `overlap`, but the genes overlapping more than 1 regions
        #   are removed.
        overlap_uniq = uniq
    )

    return(res)



def reg2gene(df, anno, no_dup = True, verbose = True):
    """Convert non-overlapping region-scale data into gene-scale.

    Each row/record in `df` will be copied (and updated accordingly) for
    each gene overlapping with the "region" within the record.

    Parameters
    ----------
    df : pandas.DataFrame
        Non-overlapping region-scale data, i.e., the regions should not
        overlap each other.
        It should contain at least four columns:
        - "chrom";
        - "start": 1-based and inclusive;
        - "end": 1-based and inclusive;
        - "region": region ID.
    anno : pandas.DataFrame
        Gene annotations.
        It should contain at least four columns:
        - "chrom";
        - "start": 1-based and inclusive;
        - "end": 1-based and inclusive;
        - "gene": gene name.
    no_dup : bool, default True
        Whether to remove genes overlapping more than 1 regions from final
        result.
    verbose : bool, default True
        Whether to show detailed logging information.

    Returns
    -------
    dict
        The convertted results.
    """
    df["chrom"] = df["chrom"].map(format_chrom)
    anno["chrom"] = anno["chrom"].map(format_chrom)

    res = get_overlap_genes(df, anno)
    overlap = res["overlap"]
    if no_dup:
         overlap = res["overlap_uniq"]

    s = "\n"
    if verbose:
        s += "\tsettings: no_dup = %s;\n" % str(no_dup)
        s += "\t#regions: unique = %d;\n" % res["n_region"]
        s += "\t#regions: having overlap genes = %d;\n" % res["n_region_overlap"]
        s += "\t#genes: unique = %d;\n" % res["n_gene"]
        s += "\t#genes: having overlap regions = %d;\n" % res["n_gene_overlap"]
        s += "\t#genes: overlap more than 1 regions = %d;\n" % res["n_gene_dup"]
        s += "\t#genes: overlap exactly 1 regions = %d;\n" % res["overlap_uniq"].shape[0]
        s += "\t#records: input = %d;\n" % df.shape[0]


    df = df.merge(overlap, on = "region", how = "inner")

    if verbose:
        s += "\t#records: output = %d;\n" % df.shape[0]
        s += "\t#regions: output unique = %d;\n" % len(df["region"].unique())
        s += "\t#genes: output unique = %d;\n" % len(df["gene"].unique())

        info("reg2gene results:\n%s" % s)


    res = dict(
        # df : pandas.DataFrame
        #   The converted data. It contains at least five columns:
        #   - "chrom";
        #   - "start": start pos of the region; 1-based and inclusive;
        #   - "end": end pos of the region; 1-based and inclusive;
        #   - "region": region ID.
        #   - "gene": name of gene overlapping the region.
        df = df,

        # overlap : pandas.DataFrame
        #   The overlapping results. It contains two columns:
        #   - "region": region ID.
        #   - "gene": name of genes overlapping the region.
        overlap = overlap
    )

    return(res)



################################################
#------------------ utils.io ------------------#
################################################

def load_cell_anno(fn, sep = "\t"):
    """Load cell annotations.

    Parameters
    ----------
    fn : str
        File storing cell annotations.
        It is a header-free file whose first two columns are:
        - "cell": cell barcode;
        - "clone": clone ID;
    sep : str, default "\t"
        Delimiter.

    Returns
    -------
    pandas.DataFrame
        Cell annotations containing two columns.
    """
    df = pd.read_csv(fn, sep = sep, header = None)
    df = df[range(2)]
    df.columns = ["cell", "clone"]
    return(df)



def load_gene_anno(fn, sep = "\t", uniq_genes = True):
    """Load gene annotations.

    Parameters
    ----------
    fn : str
        A header-free file storing gene annotations.
        Its first four columns should be:
        - "chrom";
        - "start": 1-based and inclusive;
        - "end": 1-based and inclusive;
        - "gene": gene name.
    sep : str, default "\t"
        Delimiter.
    uniq_genes : bool, default True
        Whether to only keep unique genes given their names.

    Returns
    -------
    pandas.DataFrame
        Gene annotations containing four columns.
    """
    df = pd.read_csv(fn, sep = sep, header = None)
    df = df[range(4)]
    df.columns = ["chrom", "start", "end", "gene"]
    df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].map(format_chrom)
    if uniq_genes:
        df = df.drop_duplicates("gene", ignore_index = True)
    return(df)



def load_h5ad(fn):
    return(ad.read_h5ad(fn))

def save_h5ad(adata, fn):
    return(adata.write_h5ad(fn, compression = "gzip"))



def load_truth(fn, sep = "\t"):
    """Load CNA ground truth.

    Parameters
    ----------
    fn : str
        A header-free file stroing the ground truth.
        Its first five columns should be:
        - "chrom";
        - "start": 1-based and inclusive;
        - "end": 1-based and inclusive;
        - "clone": clone ID;
        - "cna_type": CNA type, should be in {"gain", "loss", "loh"}.
    sep : str, default "\t"
        Delimiter.

    Returns
    -------
    pandas.DataFrame
        CNA ground truth.
    """
    if is_file_empty(fn):
        df = pd.DataFrame(
            columns = ["chrom", "start", "end", "clone", "cna_type"])
        for col, dtype in zip(["chrom", "start", "end", "clone", "cna_type"],
                             ["object", "int", "int", "object", "object"]):
            df[col] = df[col].astype(dtype)
        return(df)

    df = pd.read_csv(fn, sep = sep, header = None)
    df = df[range(5)]
    df.columns = ["chrom", "start", "end", "clone", "cna_type"]
    df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].map(format_chrom)
    return(df)
