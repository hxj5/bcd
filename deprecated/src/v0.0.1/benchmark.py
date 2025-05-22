# benchmark.py, converted from benchmark.R
# benchmark.R - functions for processing benchmark data.

# Note:
# 1. Several functions should be imported from `utils.R`:
#    - flush_print, str_now, write_tsv
#    - load_mtx, load_mtx3, load_gene_anno, overlap_gene_anno 
#    - reg2gene, var_reg2gene
# source utils.R
import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from datetime import datetime

# import pyreadr
# from rpy2.robjects import r
# from rpy2.robjects import pandas2ri

# pandas2ri.activate()

from utils import flush_print, str_now, write_tsv
from utils import load_mtx, load_mtx3, load_gene_anno, overlap_gene_anno 
from utils import reg2gene, var_reg2gene


def rep_n(v, n):
    func = "rep_n"
    m = len(v)
    if m < n:
        if m != 1:
            raise ValueError(f"[E::{func}] invalid length of v, 1 expected, {m} given.")
        v = v * n
    elif m != n:
        raise ValueError(f"[E::{func}] too long length of v!")
    return v

def get_step_index(step_name, step_idx):
    if step_name == "roc":
        return 5
    elif step_name == "prc":
        return 6
    else:
        raise ValueError("[E::get_step_index] invalid step name.")

def run_benchmark(sid, cnv_type, cnv_scale, dat_list, cell_anno, gene_anno, truth, out_dir,
                  overlap_mode="customize", filter_func=None, metrics=["ROC", "PRC"], max_n_cutoff=1000,
                  plot_sid=None, plot_dec=3, plot_legend_xmin=0.7, plot_legend_ymin=0.25,
                  plot_width=6.5, plot_height=5, plot_dpi=600, verbose=False, save_all=False):
    func = "run_benchmark"

    n_mt = len(metrics)
    max_n_cutoff = rep_n([max_n_cutoff], n_mt)
    if plot_sid is None:
        plot_sid = sid
    plot_dec = rep_n([plot_dec], n_mt)
    plot_legend_xmin = rep_n([plot_legend_xmin], n_mt)
    plot_legend_ymin = rep_n([plot_legend_ymin], n_mt)
    plot_width = rep_n([plot_width], n_mt)
    plot_height = rep_n([plot_height], n_mt)
    plot_dpi = rep_n([plot_dpi], n_mt)

    prefix = f"{sid}.{cnv_type}.{cnv_scale}_scale"

    # find overlap cells and genes
    flush_print(f"[I::{func}][{str_now()}] find overlap cells and genes ...")

    res_overlap = run_overlap(dat_list, mode=overlap_mode, verbose=verbose)
    cell_overlap = res_overlap['cell_overlap']
    gene_overlap = res_overlap['gene_overlap']
    dat_list = res_overlap['dat_list']

    if verbose:
        print(res_overlap)

    dir_overlap = os.path.join(out_dir, "s2_overlap")
    save_overlap(res_overlap, dir_overlap, prefix, save_all)
    flush_print(f"[I::{func}][{str_now()}] overlap data is saved to dir '{dir_overlap}'.")

    # annotate cells and genes
    flush_print(f"[I::{func}][{str_now()}] annotate cells and genes ...")

    cell_annotate = run_annotate_cell(res_overlap['cell_overlap'], cell_anno, verbose=verbose)

    if verbose:
        print(cell_annotate)

    gene_annotate = run_annotate_gene(res_overlap['gene_overlap'], gene_anno, verbose=True)

    if verbose:
        print(gene_annotate)

    dir_annotate = os.path.join(out_dir, "s3_annotate")
    save_annotate(cell_annotate, gene_annotate, dir_annotate, prefix, save_all)
    flush_print(f"[I::{func}][{str_now()}] annotate data is saved to dir '{dir_annotate}'.")

    # filter cells or genes
    if filter_func is None:
        flush_print(f"[I::{func}][{str_now()}] no filtering ...")
        cell_subset = cell_annotate['cell_anno_valid']
        gene_subset = gene_annotate['gene_anno_valid']
    else:
        flush_print(f"[I::{func}][{str_now()}] filter cells or genes ...")
        res_filter = filter_func(cell_annotate['cell_anno_valid'], gene_annotate['gene_anno_valid'])
        if verbose:
            print(res_filter)
        cell_subset = res_filter['cells']
        gene_subset = res_filter['genes']

    # subset data
    flush_print(f"[I::{func}][{str_now()}] subset data ...")
    dat_list = run_subset(cell_subset['cell'], gene_subset['Gene'], dat_list)

    if verbose:
        print(dat_list)

    dir_subset = dir_annotate
    save_subset(cell_subset, gene_subset, dat_list, dir_subset, prefix, save_all)
    flush_print(f"[I::{func}][{str_now()}] subset data is saved to dir '{dir_subset}'.")

    # construct the ground truth binary matrix
    flush_print(f"[I::{func}][{str_now()}] construct the ground truth binary matrix ...")
    res_truth = run_truth(cell_subset, gene_subset, truth, cnv_type, verbose=verbose)

    if verbose:
        print(res_truth)

    dir_truth = os.path.join(out_dir, "s4_truth")
    save_truth(res_truth, dir_truth, prefix, save_all)
    flush_print(f"[I::{func}][{str_now()}] truth data is saved to dir '{dir_truth}'.")

    # calculate metrics
    flush_print(f"[I::{func}][{str_now()}] calculate Metrics ...")

    figs = []
    s_idx = 4
    for i in range(n_mt):
        mt_name = metrics[i]
        mt_upper = mt_name.upper()
        mt_lower = mt_name.lower()
        flush_print(f"[I::{func}][{str_now()}] calculate {mt_upper} ...")

        if mt_upper == "ROC":
            dat_list = run_roc(dat_list, res_truth['truth_mtx'], max_n_cutoff=max_n_cutoff[i], strict=True, verbose=verbose)
        else:
            dat_list = run_prc(dat_list, res_truth['truth_mtx'], max_n_cutoff=max_n_cutoff[i], strict=True, verbose=verbose)

        if verbose:
            print(dat_list)

        s_idx = get_step_index(mt_lower, s_idx)
        dir_metric = os.path.join(out_dir, f"s{s_idx}_{mt_lower}")
        if mt_upper == "ROC":
            save_roc(dat_list, dir_metric, prefix, save_all)
        else:
            save_prc(dat_list, dir_metric, prefix, save_all)
        flush_print(f"[I::{func}][{str_now()}] {mt_upper} data is saved to dir '{dir_metric}'.")

        # visualization
        flush_print(f"[I::{func}][{str_now()}] visualization ...")
        p_title = None

        if cnv_type == "copy_gain":
            p_title = f"{plot_sid} {mt_upper} Curve for Copy Gain"
        elif cnv_type == "copy_loss":
            p_title = f"{plot_sid} {mt_upper} Curve for Copy Loss"
        else:
            p_title = f"{plot_sid} {mt_upper} Curve for LOH"

        if mt_upper == "ROC":
            p = plot_roc(dat_list, dec=plot_dec[i], title=p_title, legend_xmin=plot_legend_xmin[i], legend_ymin=plot_legend_ymin[i])
        else:
            p = plot_prc(dat_list, dec=plot_dec[i], title=p_title, legend_xmin=plot_legend_xmin[i], legend_ymin=plot_legend_ymin[i])

        dir_plot = dir_metric
        plot_fn = os.path.join(dir_plot, f"{prefix}.plot.{mt_lower}_figure.jpg")
        p.savefig(plot_fn, dpi=plot_dpi[i])
        flush_print(f"[I::{func}][{str_now()}] plot figure is saved to dir '{dir_plot}'.")
        figs.append(p)

    return figs

def run_bm_fast(sid, cnv_type, cnv_scale, xclone_mtx, metrics, metric_fn, cell_subset_fn, gene_subset_fn, truth_fn, out_dir,
                max_n_cutoff=1000, plot_sid=None, plot_dec=3, plot_legend_xmin=0.7, plot_legend_ymin=0.25,
                plot_width=6.5, plot_height=5, plot_dpi=600, verbose=False, save_all=False):
    func = "run_bm_fast"

    n_mt = len(metrics)
    max_n_cutoff = rep_n([max_n_cutoff], n_mt)
    if plot_sid is None:
        plot_sid = sid
    plot_dec = rep_n([plot_dec], n_mt)
    plot_legend_xmin = rep_n([plot_legend_xmin], n_mt)
    plot_legend_ymin = rep_n([plot_legend_ymin], n_mt)
    plot_width = rep_n([plot_width], n_mt)
    plot_height = rep_n([plot_height], n_mt)
    plot_dpi = rep_n([plot_dpi], n_mt)

    prefix = f"{sid}.{cnv_type}.{cnv_scale}_scale"

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # load data
    flush_print(f"[I::{func}][{str_now()}] load data ...")

    cell_subset = pd.read_csv(cell_subset_fn, sep='\t')
    gene_subset = pd.read_csv(gene_subset_fn, sep='\t')
    truth = pd.read_pickle(truth_fn)

    flush_print(f"[I::{func}][{str_now()}] the binary truth matrix is ...")
    print(truth)

    dir_input = None
    if save_all:
        dir_input = os.path.join(out_dir, "input")
        if not os.path.exists(dir_input):
            os.makedirs(dir_input)
        cell_subset.to_csv(os.path.join(dir_input, f"{prefix}.input.cell_subset.df.tsv"), sep='\t', index=False)
        gene_subset.to_csv(os.path.join(dir_input, f"{prefix}.input.gene_subset.df.tsv"), sep='\t', index=False)
        truth.to_pickle(os.path.join(dir_input, f"{prefix}.input.binary_truth.mtx.pkl"))
        flush_print(f"[I::{func}][{str_now()}] input data is saved to dir '{dir_input}'.")

    if len(cell_subset['cell']) != truth.shape[0]:
        raise ValueError(f"[E::{func}] number of subset cells and truth cells are different!")
    if not all(cell_subset['cell'].sort_values() == truth.index.sort_values()):
        raise ValueError(f"[E::{func}] subset cells and truth cells are different!")

    if len(gene_subset['Gene']) != truth.shape[1]:
        raise ValueError(f"[E::{func}] number of subset genes and truth genes are different!")
    if not all(gene_subset['Gene'].sort_values() == truth.columns.sort_values()):
        raise ValueError(f"[E::{func}] subset genes and truth genes are different!")

    # subset xclone matrix
    flush_print(f"[I::{func}][{str_now()}] subset xclone matrix ...")

    if not all(cell_subset['cell'].isin(xclone_mtx.index)):
        raise ValueError(f"[E::{func}] some subset cells are not in xclone matrix!")

    if not all(gene_subset['Gene'].isin(xclone_mtx.columns)):
        raise ValueError(f"[E::{func}] some subset genes are not in xclone matrix!")

    xclone_mtx = xclone_mtx.loc[cell_subset['cell'], gene_subset['Gene']]

    if verbose:
        print(xclone_mtx)

    # calculate metrics
    figs = []
    s_idx = 4
    for i in range(n_mt):
        mt_name = metrics[i]
        mt_upper = mt_name.upper()
        mt_lower = mt_name.lower()

        dat_list = pd.read_pickle(metric_fn[i])
        if save_all:
            dat_list.to_pickle(os.path.join(dir_input, f"{prefix}.input.{mt_lower}_list.list.pkl"))

        flush_print(f"[I::{func}][{str_now()}] calculate xclone {mt_upper} ...")

        xclone_dat_list = [{
            'sid': sid, 'cnv_type': cnv_type,
            'method': "xclone", 'method_sub': "xclone", 'mtx_type': "prob",
            'cnv_scale': cnv_scale, 'mtx': xclone_mtx, 'overlap': None
        }]

        if mt_upper == "ROC":
            xclone_dat_list = run_roc(xclone_dat_list, truth, max_n_cutoff=max_n_cutoff[i], strict=True, verbose=verbose)
        else:
            xclone_dat_list = run_prc(xclone_dat_list, truth, max_n_cutoff=max_n_cutoff[i], strict=True, verbose=verbose)

        if verbose:
            print(xclone_dat_list)

        # merge metric data
        flush_print(f"[I::{func}][{str_now()}] merge {mt_upper} data ...")

        xclone_dat = xclone_dat_list[0]

        j = 0
        for dat in dat_list:
            if dat['method'] == "xclone":
                dat_list[j] = xclone_dat
                break
            j += 1
        if j >= len(dat_list):
            dat_list.append(xclone_dat)

        if verbose:
            print(dat_list)

        s_idx = get_step_index(mt_lower, s_idx)
        dir_metric = os.path.join(out_dir, f"s{s_idx}_{mt_lower}")
        if mt_upper == "ROC":
            save_roc(dat_list, dir_metric, prefix, save_all)
        else:
            save_prc(dat_list, dir_metric, prefix, save_all)
        flush_print(f"[I::{func}][{str_now()}] updated {mt_lower} data is saved to dir '{dir_metric}'.")

        # visualization
        flush_print(f"[I::{func}][{str_now()}] visualization ...")

        if cnv_type == "copy_gain":
            p_title = f"{plot_sid} {mt_upper} Curve for Copy Gain"
        elif cnv_type == "copy_loss":
            p_title = f"{plot_sid} {mt_upper} Curve for Copy Loss"
        else:
            p_title = f"{plot_sid} {mt_upper} Curve for LOH"

        if mt_upper == "ROC":
            p = plot_roc(dat_list, dec=plot_dec[i], title=p_title, legend_xmin=plot_legend_xmin[i], legend_ymin=plot_legend_ymin[i])
        else:
            p = plot_prc(dat_list, dec=plot_dec[i], title=p_title, legend_xmin=plot_legend_xmin[i], legend_ymin=plot_legend_ymin[i])

        dir_plot = dir_metric
        plot_fn = os.path.join(dir_plot, f"{prefix}.plot.{mt_lower}_figure.jpg")
        p.savefig(plot_fn, dpi=plot_dpi[i])
        flush_print(f"[I::{func}][{str_now()}] plot figure is saved to dir '{dir_plot}'.")
        figs.append(p)

    return figs

# not modified for rds
def extract_casper(sid, dat_dir, cnv_type, method_sub="casper", mtx_type="expr", cnv_scale="gene", gene_anno=None):
    if method_sub == "casper":
        obj_fn = os.path.join(dat_dir, f"{sid}.object.pkl")
        obj = pd.read_pickle(obj_fn)
        mtx = np.log2(obj['control']['normalized']['noiseRemoved'][2])
        mtx = mtx.T  # cell x gene matrix
        dat = {'mtx': mtx, 'overlap': None}
        return dat
    elif method_sub in ["casper_median", "casper_mediandev"]:
        obj_fn = os.path.join(dat_dir, f"{sid}.object.scale.pkl")
        obj = pd.read_pickle(obj_fn)
        df = obj['segments']
        df = df[['ID', 'chr', 'start', 'end', 'state', 'median', 'medianDev']]
        df['cell'] = df['ID']
        df['chrom'] = df['chr'].str.replace("[pq]", "").str.replace("chr", "")
        if method_sub == "casper_median":
            df = df[['cell', 'chrom', 'start', 'end', 'state', 'median']].rename(columns={'median': 'score'})
            df.loc[df['state'] != "3", 'score'] = df['score'].min() - 0.1  # CHECK ME!
        else:  # method_sub == "casper_mediandev"
            df['medianDev'] = df['medianDev'].astype(float)
            df = df[['cell', 'chrom', 'start', 'end', 'state', 'medianDev']].rename(columns={'medianDev': 'score'})
            df.loc[df['state'] != "3", 'score'] = 0  # CHECK ME!
        dat = var_reg2gene(df, gene_anno)
        return dat
    else:
        raise ValueError(f"[E::extract_casper] unknown sub-type method '{method_sub}'.")


def extract_copykat(sid, dat_dir, cnv_type, method_sub="copykat", mtx_type="expr", cnv_scale="gene", gene_anno=None):
    obj_fn = os.path.join(dat_dir, f"{sid}_copykat_CNA_raw_results_gene_by_cell.txt")
    mtx = pd.read_csv(obj_fn, sep='\t', header=0, index_col=0)
    mtx = mtx.iloc[:, 7:]  # Remove the first 7 columns
    mtx = mtx.T  # Transpose to get cell x gene matrix
    dat = {'mtx': mtx, 'overlap': None}
    return dat

###############
'''
def extract_infercnv(sid, dat_dir, cnv_type, method_sub="infercnv", mtx_type="expr", cnv_scale="gene", gene_anno=None):
    obj_fn = os.path.join(dat_dir, "BayesNetOutput.HMMi6.hmm_mode-samples", "MCMC_inferCNV_obj.rds")
    obj = pd.read_pickle(obj_fn)
    mtx = obj['expr.data'].T  # cell x gene matrix
    dat = {'mtx': mtx, 'overlap': None}
    return dat


# use pyreadr to read .rds file
def extract_infercnv(sid, dat_dir, cnv_type, method_sub="infercnv", mtx_type="expr", cnv_scale="gene", gene_anno=None):
    obj_fn = os.path.join(dat_dir, "BayesNetOutput.HMMi6.hmm_mode-samples", "MCMC_inferCNV_obj.rds")
    
    # Load the RDS file using pyreadr
    try:
        result = pyreadr.read_r(obj_fn)  # Returns a dictionary where keys are object names
    except Exception as e:
        raise ValueError(f"Failed to read RDS file: {obj_fn}. Error: {e}")
    
    # Extract the 'expr.data' object from the RDS file
    if 'expr.data' in result:
        mtx = result['expr.data'].T  # Transpose the matrix (cell x gene)
    else:
        raise KeyError(f"'expr.data' not found in the RDS file: {obj_fn}")
    
    dat = {'mtx': mtx, 'overlap': None}
    return dat


# use rpy2
def extract_infercnv(sid, dat_dir, cnv_type, method_sub="infercnv", mtx_type="expr", cnv_scale="gene", gene_anno=None):
    obj_fn = os.path.join(dat_dir, "BayesNetOutput.HMMi6.hmm_mode-samples", "MCMC_inferCNV_obj.rds")
    
    # Load the RDS file using rpy2
    try:
        r_obj = r['readRDS'](obj_fn)
    except Exception as e:
        raise ValueError(f"Failed to read RDS file: {obj_fn}. Error: {e}")
    
    # Extract the expression data (assuming it's stored in 'expr.data')
    if 'expr.data' in r_obj.names:
        mtx = pandas2ri.ri2py(r_obj.rx2('expr.data')).T  # Convert to pandas DataFrame and transpose
    else:
        raise KeyError(f"'expr.data' not found in the RDS object: {obj_fn}")
    
    dat = {'mtx': mtx, 'overlap': None}
    return dat
'''

# from csv
# need to run infercnv_convert before evaluation
# sid, dat_dir, cnv_type, method_sub, mtx_type, cnv_scale, gene_anno
def extract_infercnv(sid, dat_dir, cnv_type, method_sub, mtx_type, cnv_scale, gene_anno):
    """
    Reads the matrix from a CSV file saved from the R extract_infercnv function.

    Parameters:
        sid (str): Sample ID.
        dat_dir (str): Directory where the CSV file is located.
        cnv_type (str): CNV type (e.g., "copy_loss", "copy_gain", "loh").
        output_dir (str): Directory where the CSV file is saved (optional). Defaults to dat_dir.

    Returns:
        dict: A dictionary containing the matrix and overlap (set to None).
    """

    # Construct the file path for the CSV file
    csv_file = os.path.join(dat_dir, f"{sid}_matrix.csv")
    
    # Check if the file exists
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"CSV file not found: {csv_file}")
    
    # Read the matrix from the CSV file
    try:
        mtx = pd.read_csv(csv_file, index_col=0)  # Use the first column as row indices
    except Exception as e:
        raise ValueError(f"Failed to read CSV file: {csv_file}. Error: {e}")
    
    # Prepare the dat object
    dat = {'mtx': mtx, 'overlap': None}
    return dat

def extract_numbat(sid, dat_dir, cnv_type, method_sub="numbat", mtx_type="expr", cnv_scale="gene", gene_anno=None):
    obj_fn = os.path.join(dat_dir, "joint_post_2.tsv")
    mtx = pd.read_csv(obj_fn, sep='\t')
    mtx['p_amp'] += mtx['p_bamp']
    mtx['p_del'] += mtx['p_bdel']
    # mtx['chrom'] = mtx['CHROM'].str.replace("chr", "")
    mtx['chrom'] = mtx['CHROM']
    mtx['reg_id'] = mtx.apply(lambda x: f"{x['chrom']}:{x['seg_start']}-{x['seg_end']}", axis=1)
    
    if cnv_type == "copy_gain":
        mtx = mtx.pivot(index='cell', columns='reg_id', values='p_amp')
    elif cnv_type == "copy_loss":
        mtx = mtx.pivot(index='cell', columns='reg_id', values='p_del')
    elif cnv_type == "loh":
        mtx = mtx.pivot(index='cell', columns='reg_id', values='p_loh')
    else:
        raise ValueError(f"Error: unknown cnv type '{cnv_type}'.")
    
    cells = mtx.index
    # convert to numpy array, delete to avoid issue with reg2gene
    # mtx = mtx.values
    dat = reg2gene(mtx, gene_anno)  # cell x gene matrix
    return dat

def extract_xclone(sid, dat_dir, cnv_type, method_sub="xclone", mtx_type="expr", cnv_scale="gene", gene_anno=None):
    if cnv_scale == "arm":
        mtx = load_mtx3(os.path.join(dat_dir, "cells.csv"),
                        os.path.join(dat_dir, "Features.csv"),
                        os.path.join(dat_dir, "matrix.csv"))  # cell x arm matrix
        dat = reg2gene(mtx, gene_anno)  # cell x gene matrix
    else:
        mtx = load_mtx(os.path.join(dat_dir, "cells.csv"),
                       os.path.join(dat_dir, "Features.csv"),
                       os.path.join(dat_dir, "matrix.csv"))  # cell x gene matrix
        dat = {'mtx': mtx, 'overlap': None}
    return dat

def run_extract(sid, cnv_type, cnv_scale, gene_anno, method_list, method_sub_list, mtx_type_list, dat_dir_list, verbose=False):
    func = "run_extract"
    
    if verbose:
        print(f"[I::{func}] sid = {sid}; cnv_type = {cnv_type}; cnv_scale = {cnv_scale}.")
    
    res = []
    n = len(method_list)
    for i in range(n):
        method = method_list[i].lower()
        method_sub = method_sub_list[i].lower()
        mtx_type = mtx_type_list[i].lower()
        dat_dir = dat_dir_list[i]
        
        if verbose:
            print(f"[I::{func}] processing method_sub = {method_sub}; mtx_type = {mtx_type} ...")
        
        if method == "casper":
            rdat = extract_casper(sid, dat_dir, cnv_type, method_sub, mtx_type, cnv_scale, gene_anno)
        elif method == "copykat":
            rdat = extract_copykat(sid, dat_dir, cnv_type, method_sub, mtx_type, cnv_scale, gene_anno)
        elif method == "infercnv":
            rdat = extract_infercnv(sid, dat_dir, cnv_type, method_sub, mtx_type, cnv_scale, gene_anno)
        elif method == "numbat":
            rdat = extract_numbat(sid, dat_dir, cnv_type, method_sub, mtx_type, cnv_scale, gene_anno)
        elif method == "xclone":
            rdat = extract_xclone(sid, dat_dir, cnv_type, method_sub, mtx_type, cnv_scale, gene_anno)
        else:
            raise ValueError(f"[E::{func}] unknown method '{method}'.")
        
        mtx = rdat['mtx']
        overlap = rdat['overlap']
        
        if verbose:
            print(f"[I::{func}] {mtx.shape[0]} cells and {mtx.shape[1]} genes loaded.")
        
        # remove duplicate genes
        genes = mtx.columns
        gene_cnt = genes.value_counts()
        genes1 = gene_cnt[gene_cnt == 1].index
        mtx = mtx[genes1]
        
        if verbose:
            print(f"[I::{func}] {mtx.shape[0]} cells and {mtx.shape[1]} genes left after removing duplicates.")
        
        res.append({'sid': sid, 'cnv_type': cnv_type, 'method': method, 'method_sub': method_sub, 'mtx_type': mtx_type, 'cnv_scale': cnv_scale, 'mtx': mtx, 'overlap': overlap})
    
    return res

def save_extract(dat_list, out_dir, prefix, save_all=False):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    for dat in dat_list:
        mtx = dat['mtx']
        overlap = dat['overlap']
        
        out_fn_mtx = os.path.join(out_dir, f"{dat['sid']}.{dat['cnv_type']}.{dat['method_sub']}.{dat['mtx_type']}.{dat['cnv_scale']}_scale.extract.cell_x_gene.mtx.pkl")
        mtx.to_pickle(out_fn_mtx)
        
        if overlap is not None:
            out_fn_overlap = os.path.join(out_dir, f"{dat['sid']}.{dat['cnv_type']}.{dat['method_sub']}.{dat['mtx_type']}.{dat['cnv_scale']}_scale.extract.reg2gene_mapping.df.tsv")
            overlap.to_csv(out_fn_overlap, sep='\t', index=False)
    
    dat_fn = os.path.join(out_dir, f"{prefix}.extract.data_list.list.pkl")
    pd.to_pickle(dat_list, dat_fn)

def overlap_mtx_intersect(dat_list):
    n = len(dat_list)
    for i in range(n):
        dat = dat_list[i]
        if i == 0:
            isec_cell = dat['mtx'].index
            isec_gene = dat['mtx'].columns
        else:
            isec_cell = isec_cell.intersection(dat['mtx'].index)
            isec_gene = isec_gene.intersection(dat['mtx'].columns)
    
    new_dat_list = []
    for i in range(n):
        dat = dat_list[i]
        dat['mtx'] = dat['mtx'].loc[isec_cell, isec_gene]
        dat['overlap'] = None
        new_dat_list.append(dat)
    
    return new_dat_list

def overlap_mtx_union(dat_list, na_fill=0.0):
    n = len(dat_list)
    for i in range(n):
        dat = dat_list[i]
        if i == 0:
            union_cell = dat['mtx'].index
            union_gene = dat['mtx'].columns
        else:
            union_cell = union_cell.union(dat['mtx'].index)
            union_gene = union_gene.union(dat['mtx'].columns)
    
    new_dat_list = []
    for i in range(n):
        dat = dat_list[i]
        mtx = pd.DataFrame(na_fill, index=union_cell, columns=union_gene)
        cells = dat['mtx'].index
        genes = dat['mtx'].columns
        mtx.loc[cells, genes] = dat['mtx']
        dat['mtx'] = mtx
        dat['overlap'] = None
        new_dat_list.append(dat)
    
    return new_dat_list

def overlap_mtx_customize(dat_list, na_fill=0.0):
    n = len(dat_list)
    i = 1
    j = 1
    for dat in dat_list:
        if i == 1:
            isec_cell = dat['mtx'].index
        else:
            isec_cell = isec_cell.intersection(dat['mtx'].index)
        i += 1
        if dat['method'] == "numbat":
            continue  # skip numbat as its matrix only includes target CNV regions.
        if j == 1:
            isec_gene = dat['mtx'].columns
        else:
            isec_gene = isec_gene.intersection(dat['mtx'].columns)
        j += 1
    
    new_dat_list = []
    for i in range(n):
        dat = dat_list[i]
        if dat['method'] == "numbat":
            mtx = pd.DataFrame(na_fill, index=isec_cell, columns=isec_gene)
            ic = isec_cell.intersection(dat['mtx'].index)
            ig = isec_gene.intersection(dat['mtx'].columns)
            mtx.loc[ic, ig] = dat['mtx'].loc[ic, ig]
            dat['mtx'] = mtx
        else:
            dat['mtx'] = dat['mtx'].loc[isec_cell, isec_gene]
        dat['overlap'] = None
        new_dat_list.append(dat)
    
    return new_dat_list

def run_overlap(dat_list, mode="intersect", verbose=False):
    func = "run_overlap"
    
    if mode == "intersect":
        new_dat_list = overlap_mtx_intersect(dat_list)
    elif mode == "union":
        new_dat_list = overlap_mtx_union(dat_list, na_fill=0.0)
    elif mode == "customize":
        new_dat_list = overlap_mtx_customize(dat_list, na_fill=0.0)
    else:
        raise ValueError(f"[E::{func}] unknown mode '{mode}'.")
    
    n = len(new_dat_list)
    if n <= 0:
        raise ValueError(f"[E::{func}] updated data list is empty!")
    
    for i in range(n):
        dat = new_dat_list[i]
        dat['mtx'] = dat['mtx'].fillna(0)  # na.fill?
        new_dat_list[i] = dat
    
    mtx = new_dat_list[0]['mtx']
    ovp_cells = mtx.index
    ovp_genes = mtx.columns
    
    if verbose:
        print(f"[I::{func}] shape of overlap matrix: ({len(ovp_cells)}, {len(ovp_genes)})")
    
    return {
        'cell_overlap': ovp_cells,
        'gene_overlap': ovp_genes,
        'dat_list': new_dat_list
    }

def save_overlap(res_all, out_dir, prefix, save_all=False):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    dat_list = res_all['dat_list']
    if len(dat_list) <= 0:
        return
    
    dat = dat_list[0]
    
    cell_overlap_fn = os.path.join(out_dir, f"{dat['sid']}.{dat['cnv_type']}.{dat['cnv_scale']}_scale.overlap.cells.vec.tsv")
    # Convert to a DataFrame to avoid saving index problem
    cell_overlap_df = res_all['cell_overlap'].to_frame(name='cell_overlap')
    # Save to CSV
    cell_overlap_df.to_csv(cell_overlap_fn, sep='\t', index=False, header=False)

    
    gene_overlap_fn = os.path.join(out_dir, f"{dat['sid']}.{dat['cnv_type']}.{dat['cnv_scale']}_scale.overlap.genes.vec.tsv")
    gene_overlap_df = res_all['gene_overlap'].to_frame(name='gene_overlap')
    gene_overlap_df.to_csv(gene_overlap_fn, sep='\t', index=False, header=False)
    
    if save_all:
        for dat in dat_list:
            mtx = dat['mtx']
            out_fn_mtx = os.path.join(out_dir, f"{dat['sid']}.{dat['cnv_type']}.{dat['method_sub']}.{dat['mtx_type']}.{dat['cnv_scale']}_scale.overlap.cell_x_gene.mtx.pkl")
            mtx.to_pickle(out_fn_mtx)
    
    dat_fn = os.path.join(out_dir, f"{prefix}.overlap.data_list.list.pkl")
    pd.to_pickle(dat_list, dat_fn)

def run_annotate_cell(cells, annotation, verbose=False):
    func = "run_annotate_cell"
    
    cells = pd.DataFrame({'cell': cells})
    
    cell_anno = cells.merge(annotation, on='cell', how='left')
    
    cell_anno_valid = cell_anno.dropna(subset=['cell_type'])
    cell_anno_na = cell_anno[cell_anno['cell_type'].isna()]
    
    if verbose:
        if len(cell_anno_na) > 0:
            print(f"[W::{func}] {len(cell_anno_na)} cells unannotated!")
        else:
            print(f"[I::{func}] all {len(cell_anno)} cells annotated.")
    
    return {
        'cell_anno_valid': cell_anno_valid,
        'cell_anno_na': cell_anno_na
    }

def run_annotate_gene(genes, annotation, verbose=False):
    func = "run_annotate_gene"
    
    genes = pd.DataFrame({'Gene': genes})
    
    gene_anno = genes.merge(annotation, on='Gene', how='left')
    
    gene_anno_valid = gene_anno.dropna(subset=['Chr'])
    gene_anno_na = gene_anno[gene_anno['Chr'].isna()]
    
    if verbose:
        if len(gene_anno_na) > 0:
            print(f"[W::{func}] {len(gene_anno_na)} genes unannotated!")
        else:
            print(f"[I::{func}] all {len(gene_anno)} genes annotated.")
    
    return {
        'gene_anno_valid': gene_anno_valid,
        'gene_anno_na': gene_anno_na
    }

def save_annotate(cell_annotate, gene_annotate, out_dir, prefix, save_all=False):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    cell_fn_valid = os.path.join(out_dir, f"{prefix}.annotate.valid_cells.df.tsv")
    cell_annotate['cell_anno_valid'].to_csv(cell_fn_valid, sep='\t', index=False)
    
    cell_fn_na = os.path.join(out_dir, f"{prefix}.annotate.na_cells.df.tsv")
    cell_annotate['cell_anno_na'].to_csv(cell_fn_na, sep='\t', index=False)
    
    gene_fn_valid = os.path.join(out_dir, f"{prefix}.annotate.valid_genes.df.tsv")
    gene_annotate['gene_anno_valid'].to_csv(gene_fn_valid, sep='\t', index=False)
    
    gene_fn_na = os.path.join(out_dir, f"{prefix}.annotate.na_genes.df.tsv")
    gene_annotate['gene_anno_na'].to_csv(gene_fn_na, sep='\t', index=False)

###########################################################################

def run_subset(cells, genes, dat_list):
    new_dat_list = []
    for dat in dat_list:
        dat['mtx'] = dat['mtx'].loc[cells, genes]
        new_dat_list.append(dat)
    return new_dat_list

def save_subset(cells, genes, dat_list, out_dir, prefix, save_all=False):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    if len(dat_list) <= 0:
        return
    
    dat = dat_list[0]
    
    cell_fn = os.path.join(out_dir, f"{dat['sid']}.{dat['cnv_type']}.{dat['cnv_scale']}_scale.subset.cells.df.tsv")
    cells.to_csv(cell_fn, sep='\t', index=False)
    
    gene_fn = os.path.join(out_dir, f"{dat['sid']}.{dat['cnv_type']}.{dat['cnv_scale']}_scale.subset.genes.df.tsv")
    genes.to_csv(gene_fn, sep='\t', index=False)
    
    if save_all:
        for dat in dat_list:
            mtx = dat['mtx']
            mtx_fn = os.path.join(out_dir, f"{dat['sid']}.{dat['cnv_type']}.{dat['method_sub']}.{dat['mtx_type']}.{dat['cnv_scale']}_scale.subset.cell_x_gene.mtx.pkl")
            mtx.to_pickle(mtx_fn)
    
    dat_fn = os.path.join(out_dir, f"{prefix}.subset.data_list.list.pkl")
    pd.to_pickle(dat_list, dat_fn)

def run_truth(cells, genes, truth, cnv_type, verbose=False):
    func = "run_truth"
    
    cnv_type0 = cnv_type
    truth = truth[truth['cnv_type'] == cnv_type0]
    truth['chrom'] = truth['chrom'].str.replace("chr", "")
    truth['reg_id'] = truth.apply(lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1)
    
    if verbose:
        print(f"[I::{func}] the CNV ground truth is:")
        print(truth)
    
    regions = truth[['reg_id', 'chrom', 'start', 'end']].drop_duplicates()
    
    res = overlap_gene_anno(regions, genes)
    gene_overlap = res['gene_overlap']
    n_dup = res['n_dup']
    gene_uniq = res['gene_uniq']
    
    if verbose:
        if n_dup > 0:
            print(f"[W::{func}] {n_dup} genes overlap with >1 cnv regions!")
        print(f"[I::{func}] {len(gene_uniq)} genes overlap with 1 cnv region.")
    
    truth_gene = gene_overlap.merge(truth[['reg_id', 'cell_type']], on='reg_id')
    truth_gene = truth_gene[['Gene', 'cell_type']].drop_duplicates()
    truth_gene = truth_gene.merge(cells, on='cell_type')
    truth_gene = truth_gene[['Gene', 'cell']].drop_duplicates().dropna()
    
    if verbose:
        print(f"[I::{func}] ground truth matrix mapping to gene scale:")
        print(truth_gene)
    
    truth_gene1 = truth_gene.copy()
    truth_gene1['value'] = 1
    truth_gene1 = truth_gene1.pivot(index='cell', columns='Gene', values='value').fillna(0)
    
    if verbose:
        print(f"[I::{func}] cell x gene matrix in CNV state:")
        print(truth_gene1)
    
    mtx = np.zeros((len(cells), len(genes)))
    mtx = pd.DataFrame(mtx, index=cells['cell'], columns=genes['Gene'])
    mtx.loc[truth_gene1.index, truth_gene1.columns] = truth_gene1
    
    if verbose:
        print(f"[I::{func}] dim of binary truth matrix:")
        print(mtx.shape)
    
    return {
        'gene_overlap': gene_overlap,
        'truth_mtx': mtx
    }

def save_truth(res, out_dir, prefix, save_all=False):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    gene_overlap = res['gene_overlap']
    truth_mtx = res['truth_mtx']
    
    gene_overlap_fn = os.path.join(out_dir, f"{prefix}.truth.reg2gene_mapping.df.tsv")
    gene_overlap.to_csv(gene_overlap_fn, sep='\t', index=False)
    
    mtx_fn = os.path.join(out_dir, f"{prefix}.truth.cell_x_gene.binary.mtx.pkl")
    truth_mtx.to_pickle(mtx_fn)

def run_roc(dat_list, truth_mtx, max_n_cutoff=1000, strict=True, verbose=False):
    func = "run_roc"
    if verbose:
        print(f"[I::{func}] begin ...")
    
    res = run_metric(metric="ROC", dat_list=dat_list, truth_mtx=truth_mtx, max_n_cutoff=max_n_cutoff, strict=strict, verbose=verbose)
    return res

def save_roc(dat_list, out_dir, prefix, save_all=False, verbose=False):
    func = "save_roc"
    if verbose:
        print(f"[I::{func}] begin ...")
    
    res = save_metric(metric="ROC", dat_list=dat_list, out_dir=out_dir, prefix=prefix, save_all=save_all)
    return res

def plot_roc(dat_list, dec=3, title=None, legend_xmin=0.7, legend_ymin=0.25, method_sub_case="canonical"):
    func = "plot_roc"
    res = run_plot(metric="ROC", dat_list=dat_list, dec=dec, title=title, legend_xmin=legend_xmin, legend_ymin=legend_ymin, method_sub_case=method_sub_case)
    return res

def run_prc(dat_list, truth_mtx, max_n_cutoff=1000, strict=True, verbose=False):
    func = "run_prc"
    if verbose:
        print(f"[I::{func}] begin ...")
    
    res = run_metric(metric="PRC", dat_list=dat_list, truth_mtx=truth_mtx, max_n_cutoff=max_n_cutoff, strict=strict, verbose=verbose)
    return res

def save_prc(dat_list, out_dir, prefix, save_all=False, verbose=False):
    func = "save_prc"
    if verbose:
        print(f"[I::{func}] begin ...")
    
    res = save_metric(metric="PRC", dat_list=dat_list, out_dir=out_dir, prefix=prefix, save_all=save_all)
    return res

def plot_prc(dat_list, dec=3, title=None, legend_xmin=0.7, legend_ymin=0.25, method_sub_case="canonical"):
    func = "plot_prc"
    res = run_plot(metric="PRC", dat_list=dat_list, dec=dec, title=title, legend_xmin=legend_xmin, legend_ymin=legend_ymin, method_sub_case=method_sub_case)
    return res

def run_metric(metric, dat_list, truth_mtx, max_n_cutoff=1000, strict=True, verbose=False):
    func = "run_metric"
    
    metric_upper = metric.upper()
    metric_lower = metric.lower()
    
    truth_cells = truth_mtx.index
    truth_genes = truth_mtx.columns
    
    new_dat_list = []
    
    np.random.seed(123)
    
    for dat in dat_list:
        mtx = dat['mtx']
        cells = mtx.index
        genes = mtx.columns
        dat_id = f"{dat['method_sub']}-{dat['mtx_type']}"
        
        if verbose:
            print(f"[I::{func}] begin to process {dat_id}.")
        
        if len(truth_cells) != len(cells):
            if strict:
                raise ValueError(f"[E::{func}] #cells in truth and {dat_id} matrix: {len(truth_cells)}, {len(cells)}")
            print(f"[I::{func}] #cells in truth and {dat_id} matrix: {len(truth_cells)}, {len(cells)}")
            if not all(truth_cells.isin(cells)):
                raise ValueError(f"[E::{func}] some truth-cells are not in {dat_id} matrix!")
            mtx = mtx.loc[truth_cells]
            cells = mtx.index
        
        if not all(truth_cells.sort_values() == cells.sort_values()):
            raise ValueError(f"[E::{func}] some cells in truth and {dat_id} matrix are different!")
        
        if len(truth_genes) != len(genes):
            if strict:
                raise ValueError(f"[E::{func}] #genes in truth and {dat_id} matrix: {len(truth_genes)}, {len(genes)}")
            print(f"[I::{func}] #genes in truth and {dat_id} matrix: {len(truth_genes)}, {len(genes)}")
            if not all(truth_genes.isin(genes)):
                raise ValueError(f"[E::{func}] some truth-genes are not in {dat_id} matrix!")
            mtx = mtx[truth_genes]
            genes = mtx.columns
        
        if not all(truth_genes.sort_values() == genes.sort_values()):
            raise ValueError(f"[E::{func}] some genes in truth and {dat_id} matrix are different!")
        
        if dat['mtx_type'] == "expr" and dat['cnv_type'] == "copy_loss":
            mtx = mtx * (-1)
        
        # keep the order of cell & gene the same with ground truth matrix
        mtx = mtx.loc[truth_cells, truth_genes]
        
        if verbose:
            print(f"[I::{func}] dim of final {dat_id} matrix:")
            print(mtx.shape)
        
        cutoff = np.sort(np.unique(mtx.values.flatten()))
        if len(cutoff) > max_n_cutoff:
            cutoff = np.random.choice(cutoff, size=max_n_cutoff, replace=False)
        cutoff = np.sort(np.unique(np.concatenate([cutoff, [0, 1]])))
        
        if metric_upper == "ROC":
            fpr, tpr, _ = roc_curve(truth_mtx.values.flatten(), mtx.values.flatten())
            mt_obj = {'fpr': fpr, 'tpr': tpr, 'AUC': auc(fpr, tpr)}
        else:
            precision, recall, _ = precision_recall_curve(truth_mtx.values.flatten(), mtx.values.flatten())
            mt_obj = {'precision': precision, 'recall': recall, 'AUC': auc(recall, precision)}
        
        if verbose:
            print(f"[I::{func}] AUC = {mt_obj['AUC']}.")
        
        dat['mtx'] = mtx
        dat[metric_lower] = mt_obj
        dat['auc'] = mt_obj['AUC']
        new_dat_list.append(dat)
    
    return new_dat_list

def save_metric(metric, dat_list, out_dir, prefix, save_all=False):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    metric_upper = metric.upper()
    metric_lower = metric.lower()
    
    auc_df = pd.DataFrame({
        'method': [dat['method'] for dat in dat_list],
        'method_sub': [dat['method_sub'] for dat in dat_list],
        'mtx_type': [dat['mtx_type'] for dat in dat_list],
        'auc': [dat['auc'] for dat in dat_list]
    })
    auc_fn = os.path.join(out_dir, f"{prefix}.{metric_lower}.auc.df.tsv")
    auc_df.to_csv(auc_fn, sep='\t', index=False)
    
    for dat in dat_list:
        mtx_fn = os.path.join(out_dir, f"{dat['sid']}.{dat['cnv_type']}.{dat['method_sub']}.{dat['mtx_type']}.{dat['cnv_scale']}_scale.{metric_lower}.cell_x_gene.mtx.pkl")
        if save_all and dat['mtx'] is not None:
            dat['mtx'].to_pickle(mtx_fn)
        
        metric_fn = os.path.join(out_dir, f"{dat['sid']}.{dat['cnv_type']}.{dat['method_sub']}.{dat['mtx_type']}.{dat['cnv_scale']}_scale.{metric_lower}.cardelino_cutoff.{metric_lower}.pkl")
        pd.to_pickle(dat[metric_lower], metric_fn)
    
    dat_fn = os.path.join(out_dir, f"{prefix}.{metric_lower}.pre_plot_dat_list.list.pkl")
    pd.to_pickle(dat_list, dat_fn)

def format_method_sub(method_sub, method_sub_case):
    func = "format_method_sub"
    
    if method_sub_case == "canonical":
        name_dict = {
            'casper': "CaSpER",
            'casper_median': "CaSpER_median",
            'casper_mediandev': "CaSpER_medianDev",
            'copykat': "CopyKAT",
            'infercnv': "InferCNV",
            'numbat': "Numbat",
            'xclone': "XClone"
        }
        s_str = method_sub.lower()
        if s_str not in name_dict:
            raise ValueError(f"[E::{func}] unknown method_sub '{method_sub}'.")
        t_str = name_dict[s_str]
        return t_str
    elif method_sub_case == "lower":
        name_dict = {
            'casper': "casper",
            'casper_median': "casper_median",
            'casper_mediandev': "casper_medianDev",
            'copykat': "copykat",
            'infercnv': "infercnv",
            'numbat': "numbat",
            'xclone': "xclone"
        }
        s_str = method_sub.lower()
        if s_str not in name_dict:
            raise ValueError(f"[E::{func}] unknown method_sub '{method_sub}'.")
        t_str = name_dict[s_str]
        return t_str
    elif method_sub_case == "raw":
        return method_sub
    else:
        raise ValueError(f"[E::{func}] unknown case type '{method_sub_case}'.")
    
############################################################

def format_method_sub(method_sub, method_sub_case):
    func = "format_method_sub"
    
    if method_sub_case == "canonical":
        name_dict = {
            'casper': "CaSpER",
            'casper_median': "CaSpER_median",
            'casper_mediandev': "CaSpER_medianDev",
            'copykat': "CopyKAT",
            'infercnv': "InferCNV",
            'numbat': "Numbat",
            'xclone': "XClone"
        }
        s_str = method_sub.lower()
        if s_str not in name_dict:
            raise ValueError(f"[E::{func}] unknown method_sub '{method_sub}'.")
        t_str = name_dict[s_str]
        return t_str
    elif method_sub_case == "lower":
        name_dict = {
            'casper': "casper",
            'casper_median': "casper_median",
            'casper_mediandev': "casper_medianDev",
            'copykat': "copykat",
            'infercnv': "infercnv",
            'numbat': "numbat",
            'xclone': "xclone"
        }
        s_str = method_sub.lower()
        if s_str not in name_dict:
            raise ValueError(f"[E::{func}] unknown method_sub '{method_sub}'.")
        t_str = name_dict[s_str]
        return t_str
    elif method_sub_case == "raw":
        return method_sub
    else:
        raise ValueError(f"[E::{func}] unknown case type '{method_sub_case}'.")

def run_plot(metric, dat_list, dec=3, title=None, legend_xmin=0.7, legend_ymin=0.25, method_sub_case="canonical"):
    func = "run_plot"
    
    metric_upper = metric.upper()
    metric_lower = metric.lower()
    
    p_data = pd.DataFrame()
    for dat in dat_list:
        method_sub = format_method_sub(dat['method_sub'], method_sub_case)
        mt_obj = dat[metric_lower]
        # d = pd.DataFrame(mt_obj['df'])

        if metric_lower == 'roc':
            # Extract fpr and tpr and convert to a DataFrame
            if "fpr" in mt_obj and "tpr" in mt_obj:
                df = pd.DataFrame({
                "FPR": mt_obj["fpr"],
                "TPR": mt_obj["tpr"]
                })
            else:
                raise KeyError("Either 'fpr' or 'tpr' is missing in mt_obj.")
            
        else: # prc
            if "precision" in mt_obj and "recall" in mt_obj:
                df = pd.DataFrame({
                "Precision": mt_obj["precision"],
                "Recall": mt_obj["recall"]
                })
            else:
                raise KeyError("Either 'precision' or 'recall' is missing in mt_obj.")
            
        if dec == 4:
            df['method'] = f"{method_sub}: AUC={mt_obj['AUC']:.4f}"
        else:
            df['method'] = f"{method_sub}: AUC={mt_obj['AUC']:.3f}"
        p_data = pd.concat([p_data, df], ignore_index=True)
    
    plt.figure(figsize=(6.5, 5))
    if metric_upper == "ROC":
        sns.lineplot(data=p_data, x='FPR', y='TPR', hue='method', linewidth=0.3)
        plt.xlabel("False Positive Rate (1 - Specificity)")
        plt.ylabel("True Positive Rate (Sensitivity)")
    else:
        sns.lineplot(data=p_data, x='Recall', y='Precision', hue='method', linewidth=0.3)
        plt.xlabel("Recall")
        plt.ylabel("Precision")
    
    plt.title(title)
    plt.legend(loc='lower left', bbox_to_anchor=(legend_xmin, legend_ymin), fontsize=5)
    plt.grid(False)
    plt.tight_layout()
    
    return plt