# io.py - input and output.


import anndata as ad
import pandas as pd
from .base import is_file_empty
from .grange import format_chrom



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
        return(df)

    df = pd.read_csv(fn, sep = sep, header = None)
    df = df[range(5)]
    df.columns = ["chrom", "start", "end", "clone", "cna_type"]
    df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].map(format_chrom)
    return(df)    
