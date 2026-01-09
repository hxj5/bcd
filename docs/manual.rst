Manual
======

.. contents:: Contents
   :depth: 2
   :local:



The ``bcd`` (Benchmarking of CNA Detection) pipeline evaluates the performance
of SOTA tools in CNA profile detection, tumor vs. non-tumor classification,
and subclonal structure identification.



CNA profile detection
---------------------


.. _cna-profile-quick-usage:

Quick Usage
^^^^^^^^^^^
First, please look at section :ref:`Input <cna-profile-input>` to prepare the input data.

Then call the ``cna_profile_main()`` function to run the benchmarking pipeline.

An example is:

.. code-block:: python

    from bcd.cna_profile import cna_profile_main, InferCNV, Numbat, XClone, CopyKAT, CalicoST

    infercnv = InferCNV(obj_fn = "./infercnv/BayesNetOutput.HMMi6.leiden.hmm_mode-subclusters/MCMC_inferCNV_obj.rds")
    numbat = Numbat(joint_post_fn = "./numbat/joint_post_2.tsv")
    copykat = CopyKAT(expr_mtx_fn="copykat/{sample_name}_CNA_raw_results_gene_by_cell.txt")
    xclone = XClone(combine_fn="./xclone/data/combined_final.h5ad")
    calicost = CalicoST(cnv_fn="./calicost/cnv_genelevel.tsv", clone_fn="./calicost/clone_labels.tsv")

    ret, res = cna_profile_main(
        sid = "test",
        tool_list = [infercnv, numbat],
        out_dir = "./out",
        truth_fn = "./data/truth.tsv",
        cell_anno_fn = "./data/cell_anno.tsv",
        gene_anno_fn = "./data/gene_anno.hg38.tsv",
        cna_type_list = None,         # None means ["gain", "loss", "loh"]
        verbose = True
    )
    
    print("return code = %d" % ret)


The full parameters can be found at section :ref:`Full Parameters <cna-profile-full-parameters>`.

See :ref:`Implementation <cna-profile-implementation>` for details of the pipeline.



.. _cna-profile-full-parameters:

Full Parameters
^^^^^^^^^^^^^^^

.. code-block:: python

    cna_profile_main(
        sid,
        tool_list,
        out_dir,
        truth_fn,
        cell_anno_fn,
        gene_anno_fn,
        cna_type_list = None,
        overlap_how = "isec-cells",
        max_n_cutoff = 1000,
        verbose = True
    )

    
The details are listed below:

sid : str
    Sample ID.

tool_list : list of Tool
    A list of tool-specific :class:`~.tool.Tool` objects.

out_dir : str
    The output folder.

truth_fn : str
    A header-free file stroing the ground truth.
    Its first five columns should be:

    - "chrom";
    - "start": 1-based and inclusive;
    - "end": 1-based and inclusive;
    - "clone": clone ID;
    - "cna_type": CNA type, should be in {"gain", "loss", "loh"}.

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

cna_type_list : list of str or None, default None
    A list of CNA types.
    None means using all available CNA types, including "gain",
    "loss", and "loh".

overlap_how : {"isec-cells", isec-both"}
    How to subset the tool matrices given the overlap cells and genes.

    - "isec-cells"
        Subset tool matrix by intersected cells only.
    - "isec-both"
        Subset tool matrix by intersected cells and genes.

max_n_cutoff : int or None, default 1000
    Maximum number of cutoff values for calculating metrics.
    If None, use all unique values in tool matrix.

verbose : bool, default True
    Whether to show detailed logging information.



.. _cna-profile-input:

Input
^^^^^
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



.. _cna-profile-output:

Output
^^^^^^
The final output is available at folder ``{cna_type}/3_plot``.
It contains the ROC and PRC plots for this CNA type.



.. _cna-profile-implementation:

Implementation
^^^^^^^^^^^^^^
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



Tumor vs. non-tumor classification
----------------------------------


.. _tumor-nontumor-quick-usage:

Quick Usage
^^^^^^^^^^^
First, please look at section :ref:`Input <tumor-nontumor-input>` to prepare the input data.

Then call the ``tumor_nontumor_main()`` function to run the benchmarking pipeline.

An example is:

.. code-block:: python

    from bcd.tumor_nontumor import tumor_nontumor_main, InferCNV, Numbat

    infercnv = InferCNV(obj_fn = "./infercnv/BayesNetOutput.HMMi6.leiden.hmm_mode-subclusters/MCMC_inferCNV_obj.rds")
    numbat = Numbat(clone_post_fn = "./numbat/clone_post_2.tsv")

    ret, res = tumor_nontumor_main(
        sid = "test",
        tool_list = [infercnv, numbat],
        out_dir = "./out",
        truth_fn = "./data/truth.tsv",
        tumor_labels = "cancer",
        verbose = True
    )
    
    print("return code = %d" % ret)


The full parameters can be found at section :ref:`Full Parameters <tumor-nontumor-full-parameters>`.

See :ref:`Implementation <tumor-nontumor-implementation>` for details of the pipeline.



.. _tumor-nontumor-full-parameters:

Full Parameters
^^^^^^^^^^^^^^^

.. code-block:: python

    tumor_nontumor_main(
        sid,
        tool_list,
        out_dir,
        truth_fn,
        tumor_labels,
        overlap_how = 'isec',
        fig_dpi = 300,
        verbose = True
    )
    
The details are listed below:

sid : str
    Sample ID.

tool_list : list of Tool
    A list of tool-specific :class:`~.tool.Tool` objects.

out_dir : str
    The output folder.

truth_fn : str
    A header-free file stroing the ground truth.
    Its first two columns should be:
    
    - `barcode` and `annotation`.

tumor_labels : str or list of str
    The cell type labels for tumor cells in `truth_fn`.

overlap_how : {"isec"}
    How to subset the tool matrices given the overlap cells.
    
    - "isec"
        Subset tool matrix by intersected cells only.

fig_dpi : int, default 300
    Resolution of the plot.

verbose : bool, default True
    Whether to show detailed logging information.
    


.. _tumor-nontumor-input:

Input
^^^^^
The inputs to the pipeline include:

TO BE ADDED ...



.. _tumor-nontumor-output:

Output
^^^^^^
TO BE ADDED ...



.. _tumor-nontumor-implementation:

Implementation
^^^^^^^^^^^^^^
TO BE ADDED ...



Subclonal structure identification
----------------------------------


.. _subclonal-structure-quick-usage:

Quick Usage
^^^^^^^^^^^
First, please look at section :ref:`Input <subclonal-structure-input>` to prepare the input data.

Then call the ``subclonal_structure_main()`` function to run the benchmarking pipeline.

An example is:

.. code-block:: python

    from bcd.subclonal_structure import subclonal_structure_main, InferCNV, Numbat

    infercnv = InferCNV(obj_fn = "./infercnv/BayesNetOutput.HMMi6.leiden.hmm_mode-subclusters/MCMC_inferCNV_obj.rds")
    numbat = Numbat(clone_post_fn = "./numbat/clone_post_2.tsv")

    ret, res = subclonal_structure_main(
        sid = "test",
        tool_list = [infercnv, numbat],
        out_dir = "./out",
        truth_fn = "./data/truth.tsv",
        n_cluster = 2,
        verbose = True
    )
    
    print("return code = %d" % ret)


The full parameters can be found at section :ref:`Full Parameters <subclonal-structure-full-parameters>`.

See :ref:`Implementation <subclonal-structure-implementation>` for details of the pipeline.



.. _subclonal-structure-full-parameters:

Full Parameters
^^^^^^^^^^^^^^^

.. code-block:: python

    subclonal_structure_main(
        sid,
        tool_list,
        out_dir,
        truth_fn,
        n_cluster,
        overlap_how = 'isec',
        fig_dpi = 300,
        verbose = True
    )
    
The details are listed below:

sid : str
    Sample ID.

tool_list : list of Tool
    A list of tool-specific :class:`~.tool.Tool` objects.

out_dir : str
    The output folder.

truth_fn : str
    A header-free file stroing the ground truth.
    Its first two columns should be:
    
    - `barcode` and `annotation`.

n_cluster : int
    Number of clusters for tools that do not output clone labels.

overlap_how : {"isec"}
    How to subset the tool matrices given the overlap cells.
    
    - "isec"
        Subset tool matrix by intersected cells only.

fig_dpi : int, default 300
    Resolution of the plot.

verbose : bool, default True
    Whether to show detailed logging information.
    


.. _subclonal-structure-input:

Input
^^^^^
The inputs to the pipeline include:

TO BE ADDED ...



.. _subclonal-structure-output:

Output
^^^^^^
TO BE ADDED ...



.. _subclonal-structure-implementation:

Implementation
^^^^^^^^^^^^^^
TO BE ADDED ...
