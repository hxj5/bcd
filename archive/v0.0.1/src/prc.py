# prc.py, converted from prc.R

#' Calculate PRCs
#' @param dat_list A list of CNV baf/expr/prob data for each method, see 
#'   `run_overlap` or `run_subset` for details.
#' @param truth_mtx A cell x gene binary matrix of CNV ground truth, see
#'   `run_truth` for details.
#' @param max_n_cutoff An integer. Maximum number of sampled cutoff values.
#' @param strict A bool. Whether to raise error when cell or genes in truth
#'   and method matrices are different. 
#' @param verbose A bool. Whether to output detailed information.
#' @return A list of updated CNV baf/expr/prob data.
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import precision_recall_curve, auc
from datetime import datetime

from utils import flush_print, str_now

def run_prc(dat_list, truth_mtx, max_n_cutoff=1000, strict=True, verbose=False):
    func = "run_prc"
    
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
            flush_print(f"[I::{func}] begin to process {dat_id}.")
        
        if len(truth_cells) != len(cells):
            if strict:
                raise ValueError(f"[E::{func}] #cells in truth and {dat_id} matrix: {len(truth_cells)}, {len(cells)}")
            flush_print(f"[I::{func}] #cells in truth and {dat_id} matrix: {len(truth_cells)}, {len(cells)}")
            if not all(truth_cells.isin(cells)):
                raise ValueError(f"[E::{func}] some truth-cells are not in {dat_id} matrix!")
            mtx = mtx.loc[truth_cells]
            cells = mtx.index
        
        if not all(truth_cells.sort_values() == cells.sort_values()):
            raise ValueError(f"[E::{func}] some cells in truth and {dat_id} matrix are different!")
        
        if len(truth_genes) != len(genes):
            if strict:
                raise ValueError(f"[E::{func}] #genes in truth and {dat_id} matrix: {len(truth_genes)}, {len(genes)}")
            flush_print(f"[I::{func}] #genes in truth and {dat_id} matrix: {len(truth_genes)}, {len(genes)}")
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
            flush_print(f"[I::{func}] dim of final {dat_id} matrix:")
            flush_print(mtx.shape)
        
        cutoff = np.sort(np.unique(mtx.values.flatten()))
        if len(cutoff) > max_n_cutoff:
            cutoff = np.random.choice(cutoff, size=max_n_cutoff, replace=False)
        cutoff = np.sort(cutoff)
        
        precision, recall, _ = precision_recall_curve(truth_mtx.values.flatten(), mtx.values.flatten())
        prc_auc = auc(recall, precision)
        
        if verbose:
            flush_print(f"[I::{func}] AUC = {prc_auc}.")
        
        dat['mtx'] = mtx
        dat['prc'] = {'precision': precision, 'recall': recall, 'AUC': prc_auc}
        dat['auc'] = prc_auc
        new_dat_list.append(dat)
    
    return new_dat_list

def save_prc(dat_list, out_dir, prefix, save_all=False):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    auc_df = pd.DataFrame({
        'method': [dat['method'] for dat in dat_list],
        'method_sub': [dat['method_sub'] for dat in dat_list],
        'mtx_type': [dat['mtx_type'] for dat in dat_list],
        'auc': [dat['auc'] for dat in dat_list]
    })
    auc_fn = os.path.join(out_dir, f"{prefix}.prc.auc.df.tsv")
    auc_df.to_csv(auc_fn, sep='\t', index=False)
    
    for dat in dat_list:
        mtx_fn = os.path.join(out_dir, f"{dat['sid']}.{dat['cnv_type']}.{dat['method_sub']}.{dat['mtx_type']}.{dat['cnv_scale']}_scale.prc.cell_x_gene.mtx.pkl")
        if save_all and dat['mtx'] is not None:
            dat['mtx'].to_pickle(mtx_fn)
        
        prc_fn = os.path.join(out_dir, f"{dat['sid']}.{dat['cnv_type']}.{dat['method_sub']}.{dat['mtx_type']}.{dat['cnv_scale']}_scale.prc.cardelino_cutoff.prc.pkl")
        pd.to_pickle(dat['prc'], prc_fn)
    
    dat_fn = os.path.join(out_dir, f"{prefix}.prc.pre_plot_dat_list.list.pkl")
    pd.to_pickle(dat_list, dat_fn)

def run_plot_prc(dat_list, dec=3, title=None, legend_xmin=0.7, legend_ymin=0.25, method_sub_case="canonical"):
    func = "run_plot_prc"
    
    p_data = pd.DataFrame()
    for dat in dat_list:
        method_sub = dat['method_sub']
        prc = dat['prc']
        d = pd.DataFrame({'Recall': prc['recall'], 'Precision': prc['precision']})
        if dec == 4:
            d['method'] = f"{method_sub}: AUC={prc['AUC']:.4f}"
        else:
            d['method'] = f"{method_sub}: AUC={prc['AUC']:.3f}"
        p_data = pd.concat([p_data, d])
    
    plt.figure(figsize=(6.5, 5))
    sns.lineplot(data=p_data, x='Recall', y='Precision', hue='method', size=0.3)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(title)
    plt.legend(loc='lower left', bbox_to_anchor=(legend_xmin, legend_ymin), fontsize=5)
    plt.grid(False)
    plt.tight_layout()
    
    return plt

# Main Script
work_dir = "supp_prc"
prefix = "BCH869.loh.gene_scale"
dat_list = pd.read_pickle("BCH869.loh.gene_scale.subset.data_list.list.pkl")
truth_mtx = pd.read_pickle("BCH869.loh.gene_scale.truth.cell_x_gene.binary.mtx.pkl")
sid = "BCH869"

# Core part
os.chdir(work_dir)
# Assuming benchmark.py and utils.py are already imported

func = "main_prc"

flush_print(f"[I::{func}][{str_now()}] calculate PRC ...")
dat_list = run_prc(dat_list, truth_mtx, max_n_cutoff=2000, strict=True, verbose=True)
print(dat_list)

save_prc(dat_list, work_dir, prefix, save_all=True)
flush_print(f"[I::{func}][{str_now()}] PRC data is saved to dir '{work_dir}'.")

flush_print(f"[I::{func}][{str_now()}] visualization ...")
p_title = f"{sid} PRC Curve for LOH"
p = run_plot_prc(dat_list, dec=3, title=p_title, legend_xmin=0.7, legend_ymin=0.25)

plot_fn = os.path.join(work_dir, f"{prefix}.plot.prc_figure.jpg")
p.savefig(plot_fn, dpi=600)
flush_print(f"[I::{func}][{str_now()}] plot figure is saved to dir '{work_dir}'.")

print("All Done!")