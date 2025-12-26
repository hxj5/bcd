# io.py - input and output.


import anndata as ad


def load_h5ad(fn):
    return(ad.read_h5ad(fn))

def save_h5ad(adata, fn):
    return(adata.write_h5ad(fn, compression = "gzip"))

