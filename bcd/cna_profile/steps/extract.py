# extract.py - extract CNA expression or probability matrix into adata object.


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
import scipy as sp
from logging import info
from ..utils.base import assert_e, exe_cmdline
from ..utils.io import load_gene_anno, save_h5ad



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
        info("extract matrix for '%s' ..." % tid)

        res_dir = os.path.join(out_dir, tid)
        os.makedirs(res_dir, exist_ok = True)
        
        if tid == "calicost":
            out_fn_list = [os.path.join(out_dir, "%s.%s.%s.h5ad" % \
                (out_prefix, tid, cna_type)) for cna_type in cna_type_list]
            tool.extract(
                out_fn_list = out_fn_list,
                cna_type_list = cna_type_list,
                tmp_dir = res_dir,
                verbose = verbose
            )
            for cna_type, fn in zip(cna_type_list, out_fn_list):
                out_fns[cna_type].append(fn)

        elif tid == "copykat":
            fn = os.path.join(out_dir, "%s.%s.h5ad" % (out_prefix, tid))
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
        
        elif tid == "infercnv":
            fn = os.path.join(out_dir, "%s.%s.h5ad" % (out_prefix, tid))
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
                    (out_prefix, tid, cna_type)) for cna_type in cna_type_list]
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
                    (out_prefix, tid, cna_type)) for cna_type in cna_type_list]
            tool.extract(
                out_fn_list = out_fn_list,
                cna_type_list = cna_type_list,
                tmp_dir = res_dir,
                verbose = verbose
            )
            for cna_type, fn in zip(cna_type_list, out_fn_list):
                out_fns[cna_type].append(fn)
                
        elif tid == "xclonerdr":
            out_fn_list = [os.path.join(out_dir, "%s.%s.%s.h5ad" % \
                    (out_prefix, tid, cna_type)) for cna_type in cna_type_list]
            tool.extract(
                out_fn_list = out_fn_list,
                cna_type_list = cna_type_list,
                tmp_dir = res_dir,
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
