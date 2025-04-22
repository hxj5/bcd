# truth.py - format input ground truth table into cell x gene binary matrix.


import anndata as ad
import functools
import gc
import numpy as np
import os
import pandas as pd
import scipy as sp
from logging import info, error
from ..utils.base import assert_e, expand_grid
from ..utils.gscale import reg2gene
from ..utils.io import load_truth, load_gene_anno, load_cell_anno, save_h5ad



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
                    clone, cna_type]) + "\n")
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
