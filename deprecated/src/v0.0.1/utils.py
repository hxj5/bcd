# utils.py, converted from utils.R

import pandas as pd
import numpy as np
import re
from datetime import datetime

def flush_print(s):
    print(s)
    import sys
    sys.stdout.flush()

def str_now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def write_tsv(df, file_name):
    df.to_csv(file_name, sep='\t', index=False, quoting=3)

def parse_region1(region):
    chrom = re.search(r"(.*?)(?=:|$)", region).group(0)
    
    start_match = re.search(r"(?<=:)(.*?)(?=-|$)", region)
    start = start_match.group(0).replace(",", "") if start_match else "-Inf"
    start = float(start) if start.replace('.', '', 1).isdigit() else -np.inf
    
    end_match = re.search(r"(?<=-)(.*)", region)
    end = end_match.group(0).replace(",", "") if end_match else "Inf"
    end = float(end) if end.replace('.', '', 1).isdigit() else np.inf
    
    return {"chrom": chrom, "start": start, "end": end}

def parse_regions(regions):
    res = [parse_region1(region) for region in regions]
    return pd.DataFrame(res)

def load_mtx(cell_fn, feature_fn, mtx_fn):
    cells = pd.read_csv(cell_fn, header=None)
    features = pd.read_csv(feature_fn)
    mtx = pd.read_csv(mtx_fn, header=None).values
    mtx = pd.DataFrame(mtx, index=cells[0], columns=features['GeneName'])
    return mtx

def load_mtx2(cell_fn, feature_fn, mtx_fn):
    cells = pd.read_csv(cell_fn, header=None)
    features = pd.read_csv(feature_fn)
    mtx = pd.read_csv(mtx_fn, header=None).values
    mtx = pd.DataFrame(mtx, index=cells[0], columns=[f"{chr}:{start}-{stop}" for chr, start, stop in zip(features['chr'], features['start'], features['stop'])])
    return mtx

def load_mtx3(cell_fn, feature_fn, mtx_fn):
    cells = pd.read_csv(cell_fn, header=None)
    features = pd.read_csv(feature_fn)
    mtx = pd.read_csv(mtx_fn, header=None).values
    mtx = pd.DataFrame(mtx, index=cells[0], columns=[f"{chr}:{start}-{stop}" for chr, start, stop in zip(features['chr'], features['start'], features['stop'])])
    return mtx

def load_gene_anno(gene_anno_fn):
    gene_anno = pd.read_csv(gene_anno_fn, sep='\t')
    gene_anno = gene_anno[['GeneName', 'chr', 'start', 'stop']].rename(columns={'GeneName': 'Gene', 'chr': 'Chr', 'stop': 'end'})
    gene_anno['Chr'] = gene_anno['Chr'].str.replace("chr", "")
    gene_anno = gene_anno.drop_duplicates(subset=['Gene'])
    return gene_anno

def overlap_gene_anno(regions, gene_anno):
    regions['chrom'] = regions['chrom'].str.replace("chr", "")
    gene_anno['Chr'] = gene_anno['Chr'].str.replace("chr", "")
    
    # gene_overlap = gene_anno.groupby('Gene').apply(lambda x: regions[(regions['chrom'] == x['Chr'].values[0]) & (regions['start'] <= x['end'].values[0]) & (regions['end'] >= x['start'].values[0])]['reg_id']).reset_index()
    gene_overlap = gene_anno.groupby('Gene').apply(
        lambda x: regions.loc[
            (regions['chrom'] == x['Chr'].values[0]) &
            (regions['start'] <= x['end'].values[0]) &
            (regions['end'] >= x['start'].values[0]),
            'reg_id'
    ]
    ).reset_index(name='reg_id')
    gene_overlap = gene_overlap[['Gene', 'reg_id']]
    
    gene_stat = gene_overlap.groupby('Gene').size().reset_index(name='n')
    gene_dup = gene_stat[gene_stat['n'] > 1]
    gene_uniq = gene_stat[gene_stat['n'] == 1].merge(gene_overlap, on='Gene')[['Gene', 'reg_id']]
    
    return {
        "gene_overlap": gene_overlap,
        "n_dup": len(gene_dup),
        "gene_uniq": gene_uniq
    }

def reg2gene(mtx, gene_anno, verbose=True):
    func = "reg2gene"
    
    regions = parse_regions(mtx.columns)
    regions['reg_id'] = mtx.columns
    regions = regions[['reg_id', 'chrom', 'start', 'end']]
    
    res = overlap_gene_anno(regions, gene_anno)
    gene_overlap = res['gene_overlap']
    n_dup = res['n_dup']
    gene_uniq = res['gene_uniq']
    
    if verbose:
        print(f"[I::{func}] gene_overlap:")
        print(gene_overlap)
        if n_dup > 0:
            print(f"[W::{func}] there are {n_dup} genes overlap with >1 regions!")
        print(f"[I::{func}] {len(gene_uniq)} genes overlap with 1 region.")
    
    mtx_gene = mtx[gene_uniq['reg_id']]
    mtx_gene.columns = gene_uniq['Gene']
    return {"mtx": mtx_gene, "overlap": gene_overlap}

def var_reg2gene(df, gene_anno, verbose=True):
    func = "var_reg2gene"
    
    df['reg_id'] = df.apply(lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1)
    
    gene_overlap = gene_anno.groupby('Gene').apply(lambda x: df[(df['chrom'] == x['Chr'].values[0]) & (df['start'] <= x['end'].values[0]) & (df['end'] >= x['start'].values[0])][['cell', 'reg_id', 'score']]).reset_index()
    gene_overlap.columns = ['Gene', 'cell', 'reg_id', 'score']
    
    if verbose:
        print(f"[I::{func}] gene_overlap:")
        print(gene_overlap)
    
    gene_stat = gene_overlap.groupby(['Gene', 'cell']).size().reset_index(name='n')
    gene_stat2 = gene_stat[gene_stat['n'] == 1].groupby('Gene').size().reset_index(name='m')
    
    n_cell = df['cell'].nunique()
    gene_uniq = gene_stat2[gene_stat2['m'] == n_cell].merge(gene_overlap, on='Gene')[['cell', 'Gene', 'score']]
    
    mtx_gene = gene_uniq.pivot(index='cell', columns='Gene', values='score')
    return {"mtx": mtx_gene, "overlap": gene_overlap}