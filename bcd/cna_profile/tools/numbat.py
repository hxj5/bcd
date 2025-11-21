# numbat.py


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
import scipy as sp
from logging import info
from .base import Tool
from ..utils.base import assert_e, exe_cmdline
from ..utils.gscale import reg2gene
from ..utils.io import load_gene_anno, save_h5ad


        
class Numbat(Tool):
    def __init__(self, joint_post_fn, mtx_how = 'expand'):
        """
        joint_post_fn : str
            File storing the Numbat final results.
            Typically using the "joint_post_2.tsv".
        mtx_how : {"expand", "raw"}
            How to process the extracted Numbat matrix before overlap step.
            - "expand": 
                expand the Numbat matrix to transcriptomics scale and fill value 0;
            - "raw":
                use the raw Numbat matrix.
        """
        super().__init__(
            tid = "Numbat",
            has_gain = True,
            has_loss = True,
            has_loh = True
        )
        self.joint_post_fn = joint_post_fn
        self.mtx_how = mtx_how

        
    def extract(
        self, 
        out_fn_list, 
        cna_type_list, 
        gene_anno_fn, 
        tmp_dir,
        verbose = False
    ):
        """Extract Numbat probability matrices and convert them to python objects.

        Parameters
        ----------
        out_fn_list : list of str
            Output ".h5ad" files storing the cell x gene matrix, each per cna type.
        cna_type_list : list of str
            A list of CNA types, each in {"gain", "loss", "loh"}.
        gene_anno_fn : str
            File storing gene annotations.
        tmp_dir : str
            The folder to store temporary data.
        verbose : bool, default False
            Whether to show detailed logging information.

        Returns
        -------
        Void.
        """
        return extract_cna_prob(
            joint_post_fn = self.joint_post_fn,
            out_fn_list = out_fn_list, 
            cna_type_list = cna_type_list, 
            gene_anno_fn = gene_anno_fn, 
            tmp_dir = tmp_dir,
            mtx_how = self.mtx_how,
            verbose = verbose
        )



def extract_cna_prob(
    joint_post_fn,
    out_fn_list, 
    cna_type_list, 
    gene_anno_fn,
    tmp_dir,
    mtx_how = 'expand',
    verbose = False
):
    # check args.
    if verbose:
        info("check args ...")

    assert_e(joint_post_fn)

    assert len(out_fn_list) > 0

    assert len(cna_type_list) == len(out_fn_list)
    for cna_type in cna_type_list:
        assert cna_type in ("gain", "loss", "loh")

    assert_e(gene_anno_fn)

    os.makedirs(tmp_dir, exist_ok = True)


    # load Numbat result.
    if verbose:
        info("load Numbat result ...")

    df = pd.read_csv(joint_post_fn, sep = '\t')
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
        
    return(out_fn_list)
