# extract.py - extract CNA expression or probability matrix into adata object.


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
import scipy as sp
from logging import info
from ..utils.base import assert_e, exe_cmdline
from ..utils.gscale import reg2gene
from ..utils.io import load_gene_anno, save_h5ad



def run_extract(
    args_list, out_dir, out_prefix, 
    cna_type_list, 
    gene_anno_fn, 
    numbat_mtx_how = "expand", 
    verbose = True
):
    """Extract CNA expression or probability matrix into adata object.
    
    Parameters
    ----------    
    args_list : list of ToolArgs
        A list of tool-specific :class:`~.args.ToolArgs` objects.
    out_dir : str
        The output folder.
    out_prefix : str
        The prefix to output files.
    cna_type_list : list of str
        A list of CNA types, each in {"gain", "loss", "loh"}.
    gene_anno_fn : str
        File storing gene annotations.
    numbat_mtx_how : {"expand", "raw"}
        How to process the extracted Numbat matrix before overlap step.
        - "expand": 
            expand the Numbat matrix to transcriptomics scale and fill value 0;
        - "raw":
            use the raw Numbat matrix.
    verbose : bool, default True
        Whether to show detailed logging information.
        
    Returns
    -------        
    dict
        The returned data and parameters to be used by downstream analysis.
    """
    # check args.
    assert len(args_list) > 0
    os.makedirs(out_dir, exist_ok = True)
    assert len(cna_type_list) > 0
    assert_e(gene_anno_fn)
    
    
    out_fns = {cna_type:[] for cna_type in cna_type_list}
    for args in args_list:
        tid = args.tid.lower()
        info("extract matrix for '%s' ..." % tid)

        res_dir = os.path.join(out_dir, tid)
        os.makedirs(res_dir, exist_ok = True)
        
        if tid == "infercnv":
            fn = os.path.join(out_dir, "%s.%s.h5ad" % (out_prefix, tid))
            extract_infercnv(
                obj_fn = args.obj_fn,
                out_fn = fn,
                tmp_dir = res_dir,
                verbose = verbose
            )
            for cna_type in cna_type_list:
                if args.has_cna_type(cna_type):
                    out_fns[cna_type].append(fn)
                else:
                    out_fns[cna_type].append(None)
        elif tid == "numbat":
            out_fn_list = [os.path.join(out_dir, "%s.%s.%s.h5ad" % \
                    (out_prefix, tid, cna_type)) for cna_type in cna_type_list]
            extract_numbat(
                obj_fn = args.obj_fn,
                out_fn_list = out_fn_list,
                cna_type_list = cna_type_list,
                gene_anno_fn = gene_anno_fn,
                tmp_dir = res_dir,
                mtx_how = numbat_mtx_how,
                verbose = verbose     
            )
            for cna_type, fn in zip(cna_type_list, out_fn_list):
                out_fns[cna_type].append(fn)

        elif tid == "copykat":
            out_fn_list = [os.path.join(out_dir, f"{out_prefix}.{tid}.{cna_type}.h5ad") for cna_type in cna_type_list]
            for cna_type, out_fn in zip(cna_type_list, out_fn_list):
                if args.has_cna_type(cna_type):
                    extract_copykat(
                        obj_fn=args.obj_fn,
                        out_fn=out_fn,
                        tmp_dir=res_dir,
                        verbose=verbose
                    )
                    out_fns[cna_type].append(out_fn)
                else:
                    out_fns[cna_type].append(None)

        elif tid[:6] == "xclone" and tid != "xclonerdr":
            out_fn_list = [os.path.join(out_dir, f"{out_prefix}.{tid}.{cna_type}.h5ad") for cna_type in cna_type_list]
            extract_xclone(
                obj_fn=args.obj_fn,
                out_fn_list=out_fn_list,
                cna_type_list=cna_type_list,
                tmp_dir=res_dir,
                verbose=verbose
            )
            for cna_type, fn in zip(cna_type_list, out_fn_list):
                out_fns[cna_type].append(fn)
        elif tid == "xclonerdr":
            out_fn_list = [os.path.join(out_dir, f"{out_prefix}.{tid}.{cna_type}.h5ad") for cna_type in cna_type_list]
            extract_xclone_rdr(
                obj_fn=args.obj_fn,
                rdr_fn=args.rdr_fn,
                out_fn_list=out_fn_list,
                cna_type_list=cna_type_list,
                tmp_dir=res_dir,
                verbose=verbose
            )
            for cna_type, fn in zip(cna_type_list, out_fn_list):
                out_fns[cna_type].append(fn)
        else:
            raise ValueError(f"Error: unknown tool id '{tid}'.")
                
    res = dict(
        # out_fns : dict of {str : list}
        #   Output CNA adata files.
        #   Each key is a CNA type, each value is a list of output adata files
        #   in the same order with `args_list`.
        #   If one tool does not support some CNA type, then the adata file is
        #   set to `None`.
        out_fns = out_fns
    )
    
    return(res)



def extract_infercnv(obj_fn, out_fn, tmp_dir, verbose = False):
    """Extract inferCNV expression matrix and convert it to python object.
    
    Parameters
    ----------
    obj_fn : str
        File storing the inferCNV object. Typically using the
        "BayesNetOutput.HMMi6.hmm_mode-samples/MCMC_inferCNV_obj.rds".
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



def extract_numbat(
    obj_fn, 
    out_fn_list, 
    cna_type_list, 
    gene_anno_fn, 
    tmp_dir, 
    mtx_how = "expand",
    verbose = False
):
    """Extract Numbat probability matrices and convert them to python objects.
    
    Parameters
    ----------
    obj_fn : str
        File storing the Numbat final results.
        Typically using the "joint_post_2.tsv".
    out_fn_list : list of str
        Output ".h5ad" files storing the cell x gene matrix, each per cna type.
    cna_type_list : list of str
        A list of CNA types, each in {"gain", "loss", "loh"}.
    gene_anno_fn : str
        File storing gene annotations.
    tmp_dir : str
        The folder to store temporary data.
    mtx_how : {"expand", "raw"}
        How to process the extracted Numbat matrix before overlap step.
        - "expand": 
            expand the Numbat matrix to transcriptomics scale and fill value 0;
        - "raw":
            use the raw Numbat matrix.
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
    
    assert len(out_fn_list) > 0
    
    assert len(cna_type_list) == len(out_fn_list)
    for cna_type in cna_type_list:
        assert cna_type in ("gain", "loss", "loh")
        
    assert_e(gene_anno_fn)

    os.makedirs(tmp_dir, exist_ok = True)
    
    
    # load Numbat result.
    if verbose:
        info("load Numbat result ...")
        
    df = pd.read_csv(obj_fn, sep = '\t')
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

def extract_copykat(
    obj_fn,
    out_fn,
    tmp_dir,
    verbose=False
):
    """
    cnv_type,
    method_sub="copykat",
    mtx_type="expr",
    cnv_scale="gene",
    gene_anno=None,

    Extract CNV data from CopyKAT output and save as AnnData object.

    Parameters
    ----------
    sid : str
        Sample ID.
    dat_dir : str
        Directory containing CopyKAT output.
    cnv_type : str
        CNA type (not used, for compatibility).
    method_sub : str, default "copykat"
        Sub-method name (for compatibility).
    mtx_type : str, default "expr"
        Matrix type (for compatibility).
    cnv_scale : str, default "gene"
        CNV scale (for compatibility).
    gene_anno : str or None
        Gene annotation file (not used).
    out_fn : str or None
        Output .h5ad file to save the matrix. If None, does not save.
    verbose : bool, default False
        Whether to print verbose output.

    Returns
    -------
    dict
        {"mtx": cell x gene matrix (numpy array), "overlap": None}
    """
    import pandas as pd
    import numpy as np
    import anndata as ad

    if verbose:
        info(f"Reading CopyKAT file: {obj_fn}")
    assert_e(obj_fn)
    
    os.makedirs(tmp_dir, exist_ok = True)

    mtx = pd.read_csv(obj_fn, sep="\t", header=0, dtype=str)
    mtx.index = mtx["hgnc_symbol"]
    mtx = mtx.iloc[:, 7:]  # Remove first 7 columns
    mtx = mtx.T  # cell x gene

    if verbose:
        info(f"CopyKAT matrix shape: {mtx.shape}")

    adata = ad.AnnData(
        X=mtx.values.astype(float),
        obs=pd.DataFrame(data = dict(cell = mtx.index)),
        var=pd.DataFrame(data = dict(gene = mtx.columns))
    )
    if out_fn is not None:
        save_h5ad(adata, out_fn)
        if verbose:
            info(f"Saved AnnData to {out_fn}")


def extract_xclone(
    obj_fn, 
    out_fn_list, 
    cna_type_list,
    tmp_dir,
    verbose = False
):
    """Extract XClone expression matrix and convert it to python object.
    
    Parameters
    ----------
    obj_fn : str
        File storing the XClone object.
    out_fn : str
        Output ".h5ad" file storing the cell x gene matrix.
    rdr_fn : str or None, default None
        File storing the read depth ratio (RDR) matrix.
        If None, do not use RDR.
    verbose : bool, default False
        Whether to show detailed logging information.
        
    Returns
    -------
    Void.
    """
    if verbose:
        info("check args ...")
    assert_e(obj_fn)
    assert len(out_fn_list) > 0
    assert len(cna_type_list) == len(out_fn_list)
    for cna_type in cna_type_list:
        assert cna_type in ("gain", "loss", "loh")
    os.makedirs(tmp_dir, exist_ok = True)
    if verbose:
        info("load XClone object ...")
    
    xclone_adata = ad.read_h5ad(obj_fn)
    if verbose:
        info("XClone adata shape = %s." % str(xclone_adata.shape))
    
    prob = xclone_adata.layers['prob1_merge']
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
            X=mtx,
            obs=obs_df,
            var=var_df
        )
        save_h5ad(adata, out_fn)
        
        if verbose:
            info("saved adata shape = %s." % str(adata.shape))
        
        del adata
        gc.collect()

def extract_xclone_rdr(
    obj_fn, 
    rdr_fn, 
    out_fn_list, 
    cna_type_list,
    tmp_dir,
    verbose = False
):
    """Extract XClone RDR matrix and convert it to python object.
    
    Parameters
    ----------
    obj_fn : str
        File storing the XClone object.
    rdr_fn : str
        File storing the read depth ratio (RDR) matrix.
        If None, do not use RDR.
    out_fn : str
        Output ".h5ad" file storing the cell x gene matrix.
    verbose : bool, default False
        Whether to show detailed logging information.
        
    Returns
    -------
    Void.
    """
    if verbose:
        info("check args ...")
        
    assert_e(obj_fn)
    assert_e(rdr_fn)
    
    assert len(out_fn_list) > 0
    assert len(cna_type_list) == len(out_fn_list)
    
    for cna_type in cna_type_list:
        assert cna_type in ("gain", "loss")
        
    os.makedirs(tmp_dir, exist_ok = True)
    
    if verbose:
        info("load XClone object ...")
        
    xclone_adata = ad.read_h5ad(obj_fn)
    
    if verbose:
        info("XClone adata shape = %s." % str(xclone_adata.shape))
        
    rdr_adata = ad.read_h5ad(rdr_fn)
    
    if verbose:
        info("RDR adata shape = %s." % str(rdr_adata.shape))
        
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
            mtx = xclone_adata.layers['posterior_mtx'][:, :, 0]  # copy_loss
        elif cna_type == "gain":
            mtx = xclone_adata.layers['posterior_mtx'][:, :, 1]  # copy_gain
        else:
            raise ValueError(f"Error: unknown cnv type '{cna_type}'.")

        adata = ad.AnnData(
            X=mtx,
            obs=obs_df,
            var=var_df
        )
        save_h5ad(adata, out_fn)
        if verbose:
            info("saved adata shape = %s." % str(adata.shape))
        del adata
        gc.collect()

