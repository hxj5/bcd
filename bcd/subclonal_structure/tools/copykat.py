# copykat.py


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
from logging import info, error
from logging import warning as warn
from .base import Tool
from ..utils.base import assert_e, exe_cmdline



class CopyKAT(Tool):
    def __init__(self, hclust_fn):
        """CopyKAT object.
        
        Parameters
        ----------
        hclust_fn : str
            Path to RDS file containing CopyKAT hclustering results.
        """
        super().__init__(tid = "CopyKAT")
        self.hclust_fn = hclust_fn

        
    def predict(self, out_fn, k, tmp_dir, verbose = False):
        """Process CopyKAT predictions of subclonal structures.

        Saves a TSV file with columns: 
        - ``barcode``, ``prediction``.
        """
        return predict_subclones_from_hclust(
            hclust_fn = self.hclust_fn,
            out_fn = out_fn,
            k = k,
            tmp_dir = tmp_dir,
            verbose = verbose
        )
        
        
        
def predict_subclones_from_hclust(
    hclust_fn,
    out_fn,
    k,
    tmp_dir,
    verbose = False
):
    # Check args and load data.
    assert_e(hclust_fn)
    
    s = ""
    s += '''# Extract the subclonal structure using cutree on hclust.\n'''
    s += '''\n'''
    s += '''obj <- readRDS("%s")\n''' % hclust_fn
    s += '''label <- cutree(tree = obj, k = %d)\n''' % k
    s += '''df <- data.frame(\n'''
    s += '''    barcode = gsub(".", "-", names(label), fixed = TRUE),\n'''
    s += '''    prediction = label - 1\n'''
    s += ''')\n'''
    s += '''write.table(\n'''
    s += '''    df,\n'''
    s += '''    file = "%s",\n''' % out_fn
    s += '''    sep = "\\t",\n'''
    s += '''    row.names = FALSE,\n'''
    s += '''    col.names = TRUE\n'''
    s += ''')\n'''
    s += '''\n'''

    script_fn = os.path.join(tmp_dir, "predict_subclones_from_hclust.R")
    with open(script_fn, "w") as fp:
        fp.write(s)


    # run the R script.
    if verbose:
        info("run the R script to predict subclones from hclust ...")

    exe_cmdline("Rscript %s" % script_fn)

    
    # Save to TSV
    info(f"Processed predictions saved to '{out_fn}'.")


    # Print summary
    df = pd.read_csv(out_fn, sep = '\t')
    labels = np.unique(df['prediction'])
    n_cells = len(df)
    info(f"Processed {n_cells} cells after filtering.")
    info("In total %d clones." % len(labels))
    
    return(out_fn)
