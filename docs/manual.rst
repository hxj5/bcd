Manual
======

.. contents:: Contents
   :depth: 2
   :local:



The ``bcd`` (Benchmarking of CNA Detection) pipeline evaluates the performance
of tools in detecting copy number alterations (CNAs) from single-cell and 
spatial transcriptomics, using the ROC and PRC metrics, 
given the input ground truth of CNA profiles.



Quick Usage
-----------
First, please look at section `Input`_ to prepare the input data.

Then call the ``bcd_main()`` function to run the benchmarking pipeline.

An example is:

.. code-block:: python

    from bcd import bcd_main, InferCNVArgs, NumbatArgs

    infercnv_args = InferCNVArgs(obj_fn = "./infercnv/BayesNetOutput.HMMi6.leiden.hmm_mode-subclusters/MCMC_inferCNV_obj.rds")
    numbat_args = NumbatArgs(obj_fn = "./numbat/joint_post_2.tsv")

    ret, res = bcd_main(
        sid = "test",
        args_list = [infercnv_args, numbat_args],
        out_dir = "./out",
        truth_fn = "./data/truth.tsv",
        cell_anno_fn = "./data/cell_anno.tsv",
        gene_anno_fn = "./data/gene_anno.hg38.tsv",
        cna_type_list = None,         # None means ["gain", "loss", "loh"]
        verbose = True
    )
    
    print("return code = %d" % ret)


The full parameters can be found at section `Full Parameters`_.

See `Implementation`_ for details of the pipeline.



Full Parameters
---------------

.. code-block:: python

    bcd_main(
        sid, args_list, out_dir, 
        truth_fn, cell_anno_fn, gene_anno_fn, 
        cna_type_list = None, 
        numbat_mtx_how = "expand",
        overlap_how = "isec-cells",
        max_n_cutoff = 1000,
        fig_width = 4.25, fig_height = 3.25, 
        fig_dpi = 300, fig_dec = 3, 
        fig_legend_xmin = 0.5, fig_legend_ymin = 0.12, 
        verbose = True
    )

    
The details are listed below:


sid : str
    Sample ID.
    
args_list : list of ToolArgs
    A list of tool-specific :class:`~.args.ToolArgs` objects.
    
out_dir : str
    The output folder.
        
truth_fn : str
    A header-free file stroing the ground truth.
    Its first five columns should be:
    
    * "chrom";
    * "start": 1-based and inclusive;
    * "end": 1-based and inclusive;
    * "clone": clone ID;
    * "cna_type": CNA type, should be in {"gain", "loss", "loh"}.
        
cell_anno_fn : str
    File storing cell annotations.
    It is a header-free file whose first two columns are:
    
    * "cell": cell barcode;
    * "clone": clone ID;
        
gene_anno_fn : str
    File storing gene annotations.
    It is a header-free file whose first four columns are:
    
    * "chrom";
    * "start": 1-based and inclusive;
    * "end": 1-based and inclusive;
    * "gene": gene name.
        
cna_type_list : list of str or None, default None
    A list of CNA types.
    None means using all available CNA types, including "gain",
    "loss", and "loh".
    
numbat_mtx_how : {"expand", "raw"}
    How to process the extracted Numbat matrix before overlap step.
    "expand": 
        expand the Numbat matrix to transcriptomics scale and fill value 0;
    "raw":
        use the raw Numbat matrix.

overlap_how : {"isec-cells", isec-both"}
    How to subset the tool matrices given the overlap cells and genes.
    "isec-cells"
        Subset tool matrix by intersected cells only.
    "isec-both"
        Subset tool matrix by intersected cells and genes.
        
max_n_cutoff : int or None, default 1000
    Maximum number of cutoff values for calculating metrics.
    If None, use all unique values in tool matrix.
        
fig_width : float, default 4.25
    Width of the plot in inch.
        
fig_height : float, default 3.25
    Height of the plot in inch.
        
fig_dpi : int, default 300
    Resolution of the plot.
        
fig_dec : {3, 4}
    Number of decimal places for AUC.
        
fig_legend_xmin : float, default 0.5
    The xmin position of legend.
        
fig_legend_ymin : float, default 0.12
    The ymin position of legend.
        
verbose : bool, default True
    Whether to show detailed logging information.



Input
-----
The inputs to the pipeline include:

* Tool-specific object files.
* Ground truth of clonal CNA profiles (TSV file).
* Cell annotations (TSV file).
* Feature annotations (TSV file).


Tool-specific object files
~~~~~~~~~~~~~~~~~~~~~~~~~~
The object file storing the CNA detection results of each tool, e.g.,

* for inferCNV, the ``MCMC_inferCNV_obj.rds`` file;
* for Numbat, the ``joint_post_2.tsv`` file.


Ground truth of clonal CNA profiles (TSV file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ground truth of the clonal CNA profiles, stored in a header-free TSV file.
Its first five columns ``chrom``, ``start``, ``end``, ``clone``, ``cna_type``,
where

chrom : str
    The chromosome name of the CNA region.

start : int
    The start genomic position of the CNA region, 1-based and inclusive.

end : int or "Inf"
    The end genomic position of the CNA region, 1-based and inclusive.
    To specify the end of the whole chromosome, you can use either the actual
    genomic position or simply ``Inf``.

clone : str
    The clone ID.

cna_type : str
    The CNA type, one of {"gain", "loss", "loh"}.

One clone-specific CNA per line.
An example is as follows:

.. code-block::

   chr8 1   Inf cancer1   gain
   chr6 1   Inf cancer2   loss
   chr8 1   Inf cancer2   gain
   chr6 1   Inf cancer3   loss
   chr8 1   Inf cancer3   gain
   chr11    1   Inf cancer3   loh


Cell annotations (TSV file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The cell annotation stored in a header-free TSV file.
Its first two columns are ``cell`` and ``clone``, where

cell : str
    Cell barcodes.

clone : str
    Clone ID.

An example is as follows:

.. code-block::

   AAAGATGGTCCGAAGA-1    immune
   AACCATGTCTCGTATT-1    immune
   AACGTTGTCTCTTGAT-1    cancer1
   AACTCAGAGCCTATGT-1    cancer2
   AAGACCTAGATGTAAC-1    cancer3
   AAGCCGCTCCTCAATT-1    cancer3


Feature annotations (TSV file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The feature annotation stored in a header-free TSV file.
Its first four columns are ``chrom``, ``start``, ``end``, ``feature``,
where

chrom : str
    Chromosome name of the feature.

start : int
    Start genomic position of the feature, 1-based and inclusive.

end : int
    End genomic position of the feature, 1-based and inclusive.

feature : str
    Feature name.

An example is as follows:

.. code-block::

   chr1       29554   31109   MIR1302-2HG
   chr1       34554   36081   FAM138A
   chr1       65419   71585   OR4F5
   chr2       38814   46870   FAM110C
   chr2       197569  202605  AC079779.1
   chr3       23757   24501   LINC01986



Output
------
The final output is available at folder ``{cna_type}/3_plot``.
It contains the ROC and PRC plots for this CNA type.



Implementation
--------------
The pipeline evaluates the performance of tools in detecting copy number 
alterations (CNAs) from single-cell and spatial transcriptomics, 
using the ROC and PRC metrics, given the input ground truth of CNA profiles.

It mainly includes five steps, each wrapped in one module.

The preprocessing part:

#. ``extract``: extract CNA expression or probability matrix of each tool and 
   convert the matrix into python (adata) object.
#. ``truth``: format input *clone x region* ground truth table into 
   *cell x gene* binary matrix, where entry 1 means existence of CNA.

The CNA-type-specific processing:

#. ``overlap``: subset the adata objects of tools and truth given their overlapping 
   cells and/or genes.
#. ``metric``: calculate ROC and PRC using CNA expression or probability values
   as scores and binary truth values as labels.
#. ``plot``: plot ROC and PRC.
