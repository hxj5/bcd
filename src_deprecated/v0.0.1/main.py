# main.py, converted from main.R
# Notes:
# 1. several functions should be imported from `benchmark.R`:
#    - run_benchmark, run_bm_fast.
#    - run_extract, save_extract.
# 2.  several functions should be imported from `utils.R`:
#    - flush_print, str_now.
#    - load_gene_anno.

#' Benchmark main function
#' @inheritParams run_benchmark
#' @param method_list A string vector. Names of methods.
#' @param method_sub_list A string vector. Names of sub-type methods.
#' @param mtx_type_list A string vector. Matrix type could be "baf", 
#'   "expr" or "prob".
#' @param dat_dir_list A string vector. Data dir of each method.
#' @param cell_anno_fn A string. Path to cell annotation file. It is a TSV
#'   file containing 2 columns (without header) cell and cell_type.
#' @param gene_anno_fn A string. Path to gene annotation file. It is a TSV
#'   file downloaded from XClone repo whose first 5 columns are 
#'   GeneName, GeneID, chr, start, stop.
#' @param truth_fn A string. Path to ground truth file. It is a TSV file
#'   containing at least 5 columns chrom, start, end, clone (or cell_type), 
#'   cnv_type.
#' @return A ggplot2 object. The ROC plot.

import pandas as pd
import numpy as np
import os
import re
from datetime import datetime

from benchmark import run_benchmark, run_bm_fast, run_extract, save_extract
from utils import flush_print, str_now, load_gene_anno


# Benchmark main function
def bm_main(sid, cnv_type, cnv_scale, method_list, method_sub_list, mtx_type_list, dat_dir_list,
            cell_anno_fn, gene_anno_fn, truth_fn, out_dir, overlap_mode="customize", filter_func=None, 
            metrics=["ROC", "PRC"], max_n_cutoff=1000, plot_sid=None, plot_dec=3, plot_legend_xmin=0.7, 
            plot_legend_ymin=0.25, plot_width=6.5, plot_height=5, plot_dpi=600, verbose=False, save_all=False):
    
    func = "bm_main"
    flush_print(f"[I::{func}][{str_now()}] start now ...")
    
    prefix = f"{sid}.{cnv_type}.{cnv_scale}_scale"
    
    # Load data
    flush_print(f"[I::{func}][{str_now()}] load input data ...")
    
    cell_anno = pd.read_csv(cell_anno_fn, sep='\t', header=None, names=["cell", "cell_type"])
    if verbose:
        print(cell_anno)
    
    gene_anno = load_gene_anno(gene_anno_fn)
    if verbose:
        print(gene_anno)
    
    truth_region = pd.read_csv(truth_fn, sep='\t')
    if any(~truth_region['cnv_type'].isin(["gain", "loss", "loh"])):
        raise ValueError(f"[E::{func}] cnv type should be 'gain', 'loss' or 'loh'.")
    
    idx = truth_region['cnv_type'].isin(["gain", "loss"])
    truth_region.loc[idx, 'cnv_type'] = "copy_" + truth_region.loc[idx, 'cnv_type']
    truth_region['chrom'] = truth_region['chrom'].str.replace("chr", "")
    
    if sid == "BCH869" or re.match("^GBM", sid):
        truth_region = truth_region.rename(columns={"clone": "cell_type"})
    
    truth = truth_region
    if verbose:
        print(truth)
    
    # Extract matrices
    flush_print(f"[I::{func}][{str_now()}] extract matrices ...")
    
    dat_list = run_extract(sid, cnv_type, cnv_scale, gene_anno, method_list, method_sub_list, mtx_type_list, dat_dir_list, verbose=verbose)
    
    if sid == "BCH869":
        for i, dat in enumerate(dat_list):
            cells = dat['mtx'].index
            cells = cells.str.replace("BT_", "BCH").str.replace("-", ".")
            dat['mtx'].index = cells
            dat_list[i] = dat
    
    if verbose:
        print(dat_list)
    
    dir_extract = os.path.join(out_dir, "s1_extract")
    save_extract(dat_list, dir_extract, prefix, save_all)
    flush_print(f"[I::{func}][{str_now()}] extracted data is saved to dir '{dir_extract}'.")
    
    # Benchmark core part
    flush_print(f"[I::{func}][{str_now()}] benchmark core part ...")
    
    res_bm = run_benchmark(sid, cnv_type, cnv_scale, dat_list, cell_anno, gene_anno, truth, out_dir, overlap_mode, filter_func, 
                           metrics, max_n_cutoff, plot_sid, plot_dec, plot_legend_xmin, plot_legend_ymin, plot_width, plot_height, plot_dpi, verbose, save_all)
    
    return res_bm

# Fast Benchmark (Scenario 1)
def bm_main_fast1(sid, cnv_type, cnv_scale, xclone_dir, dat_list_fn, cell_anno_fn, gene_anno_fn, truth_fn, out_dir,
                  overlap_mode="customize", filter_func=None, metrics=["ROC", "PRC"], max_n_cutoff=1000, plot_sid=None, 
                  plot_dec=3, plot_legend_xmin=0.7, plot_legend_ymin=0.25, plot_width=6.5, plot_height=5, plot_dpi=600, 
                  verbose=False, save_all=False):
    
    func = "bm_main_fast1"
    flush_print(f"[I::{func}][{str_now()}] start now ...")
    
    prefix = f"{sid}.{cnv_type}.{cnv_scale}_scale"
    
    # Load data
    flush_print(f"[I::{func}][{str_now()}] load input data ...")
    
    cell_anno = pd.read_csv(cell_anno_fn, sep='\t', header=None, names=["cell", "cell_type"])
    if verbose:
        print(cell_anno)
    
    gene_anno = load_gene_anno(gene_anno_fn)
    if verbose:
        print(gene_anno)
    
    truth_region = pd.read_csv(truth_fn, sep='\t')
    if any(~truth_region['cnv_type'].isin(["gain", "loss", "loh"])):
        raise ValueError(f"[E::{func}] cnv type should be 'gain', 'loss' or 'loh'.")
    
    idx = truth_region['cnv_type'].isin(["gain", "loss"])
    truth_region.loc[idx, 'cnv_type'] = "copy_" + truth_region.loc[idx, 'cnv_type']
    truth_region['chrom'] = truth_region['chrom'].str.replace("chr", "")
    
    if sid == "BCH869" or re.match("^GBM", sid):
        truth_region = truth_region.rename(columns={"clone": "cell_type"})
    
    truth = truth_region
    if verbose:
        print(truth)
    
    # Extract matrices
    flush_print(f"[I::{func}][{str_now()}] extract XClone matrices ...")
    
    xclone_dat_list = run_extract(sid, cnv_type, cnv_scale, gene_anno, ["xclone"], ["xclone"], ["prob"], [xclone_dir], verbose)
    xclone_dat = xclone_dat_list[0]
    
    if sid == "BCH869":
        cells = xclone_dat['mtx'].index
        cells = cells.str.replace("BT_", "BCH").str.replace("-", ".")
        xclone_dat['mtx'].index = cells
    
    if verbose:
        print(xclone_dat)
    
    # Merge data
    flush_print(f"[I::{func}][{str_now()}] merge matrices data ...")
    
    dat_list = pd.read_pickle(dat_list_fn)
    for i, dat in enumerate(dat_list):
        if dat['method'] == "xclone":
            dat_list[i] = xclone_dat
            break
    else:
        dat_list.append(xclone_dat)
    
    if verbose:
        print(dat_list)
    
    dir_extract = os.path.join(out_dir, "s1_extract")
    save_extract(dat_list, dir_extract, prefix, save_all)
    flush_print(f"[I::{func}][{str_now()}] extracted data is saved to dir '{dir_extract}'.")
    
    # Benchmark core part
    flush_print(f"[I::{func}][{str_now()}] benchmark core part ...")
    
    res_bm = run_benchmark(sid, cnv_type, cnv_scale, dat_list, cell_anno, gene_anno, truth, out_dir, overlap_mode, filter_func, 
                           metrics, max_n_cutoff, plot_sid, plot_dec, plot_legend_xmin, plot_legend_ymin, plot_width, plot_height, plot_dpi, verbose, save_all)
    
    return res_bm

# Fast Benchmark (Scenario 2)
def bm_main_fast2(sid, cnv_type, cnv_scale, xclone_dir, metrics, metric_fn, gene_anno_fn, cell_subset_fn, gene_subset_fn, truth_fn, out_dir,
                  overlap_mode="customize", filter_func=None, max_n_cutoff=1000, plot_sid=None, plot_dec=3, plot_legend_xmin=0.7, 
                  plot_legend_ymin=0.25, plot_width=6.5, plot_height=5, plot_dpi=600, verbose=False, save_all=False):
    
    func = "bm_main_fast2"
    flush_print(f"[I::{func}][{str_now()}] start now ...")
    
    prefix = f"{sid}.{cnv_type}.{cnv_scale}_scale"
    
    # Load data
    flush_print(f"[I::{func}][{str_now()}] load input data ...")
    
    gene_anno = load_gene_anno(gene_anno_fn)
    if verbose:
        print(gene_anno)
    
    # Extract matrices
    flush_print(f"[I::{func}][{str_now()}] extract XClone matrices ...")
    
    xclone_dat_list = run_extract(sid, cnv_type, cnv_scale, gene_anno, ["xclone"], ["xclone"], ["prob"], [xclone_dir], verbose)
    xclone_dat = xclone_dat_list[0]
    xclone_mtx = xclone_dat['mtx']
    
    if sid == "BCH869":
        cells = xclone_mtx.index
        cells = cells.str.replace("BT_", "BCH").str.replace("-", ".")
        xclone_mtx.index = cells
    
    if verbose:
        print(xclone_dat)
    
    # Benchmark core part
    flush_print(f"[I::{func}][{str_now()}] benchmark core part ...")
    
    res_bm = run_bm_fast(sid, cnv_type, cnv_scale, xclone_mtx, metrics, metric_fn, cell_subset_fn, gene_subset_fn, truth_fn, out_dir, 
                         max_n_cutoff, plot_sid, plot_dec, plot_legend_xmin, plot_legend_ymin, plot_width, plot_height, plot_dpi, verbose, save_all)
    
    return res_bm

def construct_truth(truth_region, cnv_cell_type):
    truth = pd.DataFrame()
    for ct in cnv_cell_type:
        truth = pd.concat([truth, truth_region.assign(cell_type=ct)], ignore_index=True)
    return truth