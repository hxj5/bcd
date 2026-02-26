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

    infercnv = InferCNV(obj_fn = "./infercnv/run.final.infercnv_obj")
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


Multiple runs of the same tool
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To compare multiple runs from the same tool (e.g., different random seeds or parameters),
pass ``run_id`` when creating each tool instance. Each run will appear as a separate
curve in ROC/PRC plots and in the legend:

.. code-block:: python

    from bcd.cna_profile import cna_profile_main, CopyKAT, InferCNV

    # Multiple CopyKAT runs with different parameters
    copykat_run1 = CopyKAT(expr_mtx_fn="copykat/run1_CNA_raw_results.txt", run_id="run1")
    copykat_run2 = CopyKAT(expr_mtx_fn="copykat/run2_CNA_raw_results.txt", run_id="run2")

    # Multiple InferCNV runs
    infercnv_rep1 = InferCNV(obj_fn="infercnv/rep1.rds", run_id="rep1")
    infercnv_rep2 = InferCNV(obj_fn="infercnv/rep2.rds", run_id="rep2")

    ret, res = cna_profile_main(
        sid="test",
        tool_list=[copykat_run1, copykat_run2, infercnv_rep1, infercnv_rep2],
        out_dir="./out",
        truth_fn="./data/truth.tsv",
        cell_anno_fn="./data/cell_anno.tsv",
        gene_anno_fn="./data/gene_anno.hg38.tsv",
    )


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
    For multiple runs of the same tool, pass ``run_id`` to each instance
    (e.g., ``CopyKAT(expr_mtx_fn="...", run_id="run1")``) to distinguish
    them in plots and output files.

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
The object file storing the CNA detection results of each tool. Below are the input requirements for each tool:

**InferCNV**
    - Input: RDS file (``obj_fn``) storing the inferCNV object (typically ``run.final.infercnv_obj``)
    - Format: R Seurat object with expression matrix (``obj@expr.data``)
    - Output: Cell × gene CNA expression matrix used for ROC/PRC analysis

**Numbat**
    - Input: TSV file (``joint_post_fn``) with genomic segments and CNV posterior probabilities
    - Format: Tab-delimited file with columns: chrom, start, end, cell, cnv_status, p_cnv, posterior
    - Output: Cell × region CNA probability matrix

**CopyKAT**
    - Input: TSV/CSV file (``expr_mtx_fn``) with raw CNA results (gene × cell matrix)
    - Format: Typically named ``{sample_name}_CNA_raw_results_gene_by_cell.txt``
    - Output: Gene × cell CNA expression matrix
    - Notes: Can contain relative copy number values or binary CNA calls

**XClone**
    - Input: H5AD file (``combine_fn``) with combined CNA analysis results
    - Format: AnnData object with CNA probability data stored in ``.layers['prob1_merge_refined']``
    - Typical path: ``./xclone/data/combined_final.h5ad``
    - Output: Cell × gene CNA probability matrix

**CalicoST**
    - Input: Two TSV files
    
      - CNV file (``cnv_fn``): Gene-level CNV calls (typically ``cnv_genelevel.tsv``)
      - Clone file (``clone_fn``): Clone labels and assignments (typically ``clone_labels.tsv``)
    
    - Format: Tab-delimited files
    - CNV file columns: chrom, start, end, gene, clone, cnv_type (or similar)
    - Clone file columns: cell, clone
    - Output: Cell × gene CNA matrix

All tools support an optional ``run_id`` parameter for multiple runs
(e.g., ``CopyKAT(expr_mtx_fn="...", run_id="run1")``).


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

    from bcd.tumor_nontumor import (
        tumor_nontumor_main,
        InferCNV,
        Numbat,
        CopyKAT,
        XClone,
        CalicoST
    )

    # Define arguments for each tool
    infercnv = InferCNV(obj_fn="./infercnv/run.final.infercnv_obj")
    numbat = Numbat(clone_post_fn="./numbat/clone_post_2.tsv")
    copykat = CopyKAT(ploidy_pred_fn="./copykat/copykat_pred.tsv")
    xclone = XClone(xclone_tumor_pred_fn="./xclone/xclone_tumor_predictions.tsv")
    calicost = CalicoST(tumor_prop_fn="./calicost/tumor_proportions.tsv")

    # Run tumor vs normal classification
    ret, res = tumor_nontumor_main(
        sid="HCC-3",
        tool_list=[infercnv, numbat, copykat, xclone],
        out_dir="./out",
        truth_fn="./data/truth.tsv",
        tumor_labels="tumor",
        overlap_how="isec",
        fig_dpi=300,
        verbose=True
    )
    
    print("return code = %d" % ret)


Multiple runs of the same tool
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To compare multiple runs from the same tool, pass ``run_id`` when creating
each tool instance. Each run will appear as a separate subplot and in the
metrics bar chart:

.. code-block:: python

    from bcd.tumor_nontumor import tumor_nontumor_main, CopyKAT, InferCNV

    copykat_run1 = CopyKAT(ploidy_pred_fn="./copykat/run1_pred.tsv", run_id="run1")
    copykat_run2 = CopyKAT(ploidy_pred_fn="./copykat/run2_pred.tsv", run_id="run2")
    infercnv_rep1 = InferCNV(obj_fn="./infercnv/rep1.rds", run_id="rep1")
    infercnv_rep2 = InferCNV(obj_fn="./infercnv/rep2.rds", run_id="rep2")

    ret, res = tumor_nontumor_main(
        sid="HCC-3",
        tool_list=[copykat_run1, copykat_run2, infercnv_rep1, infercnv_rep2],
        out_dir="./out",
        truth_fn="./data/truth.tsv",
        tumor_labels="tumor",
    )


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
    For multiple runs of the same tool, pass ``run_id`` to each instance
    (e.g., ``CopyKAT(ploidy_pred_fn="...", run_id="run1")``) to distinguish
    them in plots and output files.

out_dir : str
    The output folder.

truth_fn : str
    A header-free file stroing the ground truth.
    Its first two columns should be:
    
    - `barcode` and `annotation`.

tumor_labels : str or list of str
    The cell type labels for tumor cells in `truth_fn`.
    Can be a single label (e.g., "tumor") or a list of labels
    (e.g., ["cancer", "malignant"]) that should be classified as tumor.

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

The pipeline requires the following inputs:

1. **Tool-specific files**: Predictions/outputs from each tool in the ``tool_list``.

2. **Ground truth file** (``truth_fn``): A header-free TSV file with two columns:
   
   - ``barcode``: Cell barcode/ID
   - ``annotation``: Cell type annotation (any string labels)
   
   Example:
   
   .. code-block:: text
   
       AAACCCAAGAAACTG  tumor
       AAACCCAAGAAACTG  normal
       AAACCCAAGAAACTT  tumor
       AAACCCAAGAAACTT  normal
       ...

Tool-specific input files are described below:

**InferCNV**
    - Input: RDS file (``obj_fn``) storing the inferCNV object (typically ``run.final.infercnv_obj``)
    - Format: R Seurat object with expression matrix

**Numbat**
    - Input: TSV file (``clone_post_fn``) with columns: ``cell`` and ``compartment_opt``
    - Format: Tab-delimited file
    - Example columns: cell, chrom, p_cnv, state, ..., compartment_opt (normal/tumor)

**CopyKAT**
    - Input: TSV file (``ploidy_pred_fn``) with columns: ``cell.names`` and ``copykat.pred``
    - Format: Tab-delimited or CSV
    - Valid predictions: "diploid" (normal) or "aneuploid" (tumor)

**XClone**
    - Input: TSV file (``xclone_tumor_pred_fn``) with columns: ``barcode`` and ``prediction``
    - Format: Tab-delimited with header
    - Valid predictions: "normal" or "tumor"

**CalicoST**
    - Input: TSV file (``tumor_prop_fn``) with columns: ``BARCODES``, ``clone_label``, ``tumor_proportion``
    - Format: Tab-delimited
    - Uses K-means clustering on ``tumor_proportion`` to classify cells

All tools support an optional ``run_id`` parameter for multiple runs
(e.g., ``CopyKAT(ploidy_pred_fn="...", run_id="run1")``).



.. _tumor-nontumor-output:

Output
^^^^^^

The pipeline generates outputs organized in the following directory structure:

.. code-block:: text

    out_dir/
    ├── 0_pp/                           # Preprocessing
    │   ├── <tool_id>/
    │   │   └── <tool_id>_predictions.tsv
    │   └── truth/
    │       └── truth.tsv
    ├── 1_overlap/                      # Cell overlap analysis
    │   ├── overlap.tools.intersect.cells.tsv
    │   ├── overlap.tools_and_truth.intersect.cells.tsv
    │   ├── overlap.<tool_id>.tsv       # Per-tool subset (1 per tool)
    │   └── overlap.truth.tsv
    ├── 2_metric/
    │   └── metrics.tsv                 # Evaluation metrics
    └── 3_plot/
        ├── <sample_id>.labels.confusion_matrix.jpg
        ├── <sample_id>.labels.histogram.jpg
        ├── <sample_id>.metrics.bar.jpg
        └── <sample_id>.metrics.radar.jpg

**Output File Descriptions**

``metrics.tsv``
    Performance metrics for each tool. Columns:
    
    - ``tool``: Tool name
    - ``metric``: Metric type (accuracy, precision, recall, F1, ARI)
    - ``value``: Metric value (0-1 range)

``overlap.<tool_id>.tsv``
    Subset of each tool's predictions, including only cells present in all tools and ground truth.
    Columns: ``barcode``, ``prediction``, plus tool-specific columns.

``overlap.truth.tsv``
    Ground truth labels for overlapping cells.

**Visualization Files**

- ``confusion_matrix.jpg``: Confusion matrices for all tools in a grid layout
- ``histogram.jpg``: Distribution of predicted tumor/normal labels per tool
- ``metrics.bar.jpg``: Bar plots of metrics (accuracy, precision, recall, F1, ARI) per tool
- ``metrics.radar.jpg``: Radar chart showing all metrics for each tool



.. _tumor-nontumor-implementation:

Implementation
^^^^^^^^^^^^^^

The pipeline performs the following steps:

1. **Preprocessing** (Stage 0)
   
   - Each tool is processed to extract tumor predictions
   - Outputs are standardized to TSV format with columns: ``barcode``, ``prediction``
   - Tool-specific processing:
   
     - **InferCNV**: Uses hierarchical clustering on expression matrix with 2 clusters (normal/tumor)
     - **Numbat**: Extracts ``compartment_opt`` column
     - **CopyKAT**: Maps "diploid" → "normal", "aneuploid" → "tumor"
     - **XClone**: Validates input predictions
     - **CalicoST**: Applies K-means clustering on tumor_proportion values

2. **Cell Overlap Analysis** (Stage 1)
   
   - Identifies cells present in ALL tools and ground truth
   - Generates intersection lists for quality control
   - Subsets all prediction files to contain only overlapping cells
   - Ensures fair comparison across all tools

3. **Metrics Calculation** (Stage 2)
   
   - Converts binary predictions using positive label "tumor"
   - Calculates for each tool:
   
     - **Accuracy**: (TP + TN) / (TP + TN + FP + FN)
     - **Precision**: TP / (TP + FP)
     - **Recall**: TP / (TP + FN)
     - **F1 Score**: 2 × (Precision × Recall) / (Precision + Recall)
     - **ARI** (Adjusted Rand Index): Clustering similarity metric

4. **Visualization** (Stage 3)
   
   - **Confusion matrices**: Shows prediction vs ground truth for each tool
   - **Histograms**: Distribution of predicted labels
   - **Bar plots**: Performance metrics comparison
   - **Radar chart**: Multi-metric performance visualization

**Key Features**

- **Standardized input format**: All tools output predictions to TSV with ``barcode`` and ``prediction`` columns
- **Cell overlap tracking**: Ensures valid comparisons by using common cell set
- **Comprehensive evaluation**: Multiple metrics for robust performance assessment
- **Automated visualization**: Generates publication-quality figures



Subclonal structure identification
----------------------------------


.. _subclonal-structure-quick-usage:

Quick Usage
^^^^^^^^^^^
First, please look at section :ref:`Input <subclonal-structure-input>` to prepare the input data.

Then call the ``subclonal_structure_main()`` function to run the benchmarking pipeline.

An example is:

.. code-block:: python

    from bcd.subclonal_structure import (
        subclonal_structure_main,
        CalicoST,
        CopyKAT,
        InferCNV,
        Numbat,
        XClone
    )

    # Define arguments for each tool
    calicost = CalicoST(clone_label_fn="./calicost/clone_labels.tsv")
    copykat = CopyKAT(hclust_fn="./copykat/hclust.rds")
    infercnv = InferCNV(obj_fn="./infercnv/run.final.infercnv_obj")
    numbat = Numbat(clone_post_fn="./numbat/clone_post_2.tsv")
    xclone = XClone(clone_post_fn="./xclone/HCC3_select_ref_cell_clone_posteriors.tsv")

    # Run subclonal structure identification
    ret, res = subclonal_structure_main(
        sid="test",
        tool_list=[calicost, copykat, infercnv, numbat, xclone],
        out_dir="./out",
        truth_fn="./data/truth.tsv",
        n_cluster=3,  # Number of clusters for tools that need clustering (CopyKAT, InferCNV)
        overlap_how="isec",
        fig_dpi=300,
        verbose=True
    )
    
    print("return code = %d" % ret)


Multiple runs of the same tool
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To compare multiple runs from the same tool, pass ``run_id`` when creating
each tool instance. Each run will appear as a separate subplot in confusion
matrices and in the metrics bar chart:

.. code-block:: python

    from bcd.subclonal_structure import subclonal_structure_main, CopyKAT, InferCNV

    copykat_run1 = CopyKAT(hclust_fn="./copykat/run1_hclust.rds", run_id="run1")
    copykat_run2 = CopyKAT(hclust_fn="./copykat/run2_hclust.rds", run_id="run2")
    infercnv_rep1 = InferCNV(obj_fn="./infercnv/rep1.rds", run_id="rep1")
    infercnv_rep2 = InferCNV(obj_fn="./infercnv/rep2.rds", run_id="rep2")

    ret, res = subclonal_structure_main(
        sid="test",
        tool_list=[copykat_run1, copykat_run2, infercnv_rep1, infercnv_rep2],
        out_dir="./out",
        truth_fn="./data/truth.tsv",
        n_cluster=3,
    )


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
    For multiple runs of the same tool, pass ``run_id`` to each instance
    (e.g., ``CopyKAT(hclust_fn="...", run_id="run1")``) to distinguish
    them in plots and output files.

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
The pipeline requires the following inputs:

1. **Tool-specific files**: Predictions/outputs from each tool in the ``tool_list``.

2. **Ground truth file** (``truth_fn``): A header-free TSV file with two columns:
   
   - ``barcode``: Cell barcode/ID
   - ``annotation``: Clone/subclone label (any string labels)
   
   Example:
   
   .. code-block:: text
   
       AAACCCAAGAAACTG  clone1
       AAACCCAAGAAACTG  clone2
       AAACCCAAGAAACTT  clone1
       AAACCCAAGAAACTT  clone3
       ...

Tool-specific input files are described below:

**CalicoST**
    - Input: TSV file (``clone_label_fn``) with columns: ``BARCODES`` and ``clone_label``
    - Format: Tab-delimited file
    - Output: Extracts clone labels and maps them to numeric predictions
    - Notes: Rows with empty ``clone_label`` values are filtered out

**CopyKAT**
    - Input: RDS file (``hclust_fn``) containing CopyKAT hierarchical clustering results
    - Format: R object storing hclust tree
    - Output: Uses ``cutree()`` to cut the tree into ``k`` clusters
    - Notes: Requires ``n_cluster`` parameter to specify number of clusters

**InferCNV**
    - Input: RDS file (``obj_fn``) storing the inferCNV object (typically ``run.final.infercnv_obj``)
    - Format: R inferCNV object with hierarchical clustering stored in ``obj@tumor_subclusters$hc``
    - Output: Cuts existing hierarchical clustering trees into ``k`` clusters per group
    - Notes: 
      - Requires ``n_cluster`` parameter to specify number of clusters per group
      - Must have run ``infercnv::run(analysis_mode='subclusters')`` to generate clustering
      - Processes each group separately and combines results

**Numbat**
    - Input: TSV file (``clone_post_fn``) with columns: ``cell`` and ``clone_opt``
    - Format: Tab-delimited file
    - Output: Extracts clone labels from ``clone_opt`` column
    - Notes: Rows with empty ``clone_opt`` values are filtered out

**XClone**
    - Input: TSV file (``clone_post_fn``) with columns: ``cell_barcode`` and ``clone_id_refined``
    - Format: Tab-delimited file
    - Output: Extracts clone labels from ``clone_id_refined`` column
    - Notes: 
      - Barcode column defaults to ``cell_barcode`` but can auto-detect alternatives (``cell``, ``barcode``, etc.)
      - Rows with empty ``clone_id_refined`` values are filtered out

All tools support an optional ``run_id`` parameter for multiple runs
(e.g., ``CopyKAT(hclust_fn="...", run_id="run1")``).



.. _subclonal-structure-output:

Output
^^^^^^
The pipeline generates outputs organized in the following directory structure:

.. code-block:: text

    out_dir/
    ├── 0_pp/                           # Preprocessing
    │   ├── <tool_id>/
    │   │   └── <tool_id>_predictions.tsv
    │   └── truth/
    │       └── truth.tsv
    ├── 1_overlap/                      # Cell overlap analysis
    │   ├── overlap.tools.intersect.cells.tsv
    │   ├── overlap.tools_and_truth.intersect.cells.tsv
    │   ├── overlap.<tool_id>.tsv       # Per-tool subset (1 per tool)
    │   └── overlap.truth.tsv
    ├── 2_metric/
    │   └── metrics.tsv                 # Evaluation metrics
    └── 3_plot/
        ├── <sample_id>.labels.confusion_matrix.jpg
        └── <sample_id>.metrics.bar.jpg

**Output File Descriptions**

``<tool_id>_predictions.tsv``
    Preprocessed predictions for each tool. Columns:
    
    - ``barcode``: Cell barcode
    - ``prediction``: Numeric clone ID (0, 1, 2, ...)
    - ``clone_label`` (or ``cluster_id``): Original clone label from tool

``truth.tsv``
    Formatted ground truth labels. Columns:
    
    - ``barcode``: Cell barcode
    - ``annotation``: Numeric clone ID (mapped from original labels)
    - ``annotation_old``: Original clone label

``metrics.tsv``
    Performance metrics for each tool. Columns:
    
    - ``tool``: Tool name
    - ``metric``: Metric type (currently only ``ARI``)
    - ``value``: Metric value (0-1 range for ARI)

``overlap.<tool_id>.tsv``
    Subset of each tool's predictions, including only cells present in all tools and ground truth.
    Columns: ``barcode``, ``prediction``, plus tool-specific columns (e.g., ``clone_label``).

``overlap.truth.tsv``
    Ground truth labels for overlapping cells.

**Visualization Files**

- ``labels.confusion_matrix.jpg``: Confusion matrices for all tools in a grid layout, showing predicted vs. true clone assignments
- ``metrics.bar.jpg``: Bar plots of ARI (Adjusted Rand Index) for each tool



.. _subclonal-structure-implementation:

Implementation
^^^^^^^^^^^^^^
The pipeline evaluates the performance of tools in identifying subclonal structures from single-cell and spatial transcriptomics data, using the Adjusted Rand Index (ARI) metric.

It mainly includes four steps, each wrapped in one module:

1. **Preprocessing** (Stage 0)
   
   - Each tool is processed to extract clone/subclone predictions
   - Outputs are standardized to TSV format with columns: ``barcode``, ``prediction``, and original label
   - Tool-specific processing:
   
     - **CalicoST**: Extracts ``clone_label`` column, filters empty values, maps to numeric IDs
     - **CopyKAT**: Uses ``cutree()`` on hierarchical clustering tree with ``k`` clusters
     - **InferCNV**: Cuts existing hierarchical clustering trees (from ``obj@tumor_subclusters$hc``) into ``k`` clusters per group, then combines results
     - **Numbat**: Extracts ``clone_opt`` column, filters empty values, maps to numeric IDs
     - **XClone**: Extracts ``clone_id_refined`` column, filters empty values, maps to numeric IDs
   
   - Ground truth is formatted: original labels are mapped to numeric IDs (0, 1, 2, ...)

2. **Cell Overlap Analysis** (Stage 1)
   
   - Identifies cells present in ALL tools and ground truth
   - Generates intersection lists for quality control:
     - ``overlap.tools.intersect.cells.tsv``: Cells present in all tools
     - ``overlap.tools_and_truth.intersect.cells.tsv``: Cells present in all tools AND ground truth
   - Subsets all prediction files to contain only overlapping cells
   - Ensures fair comparison across all tools

3. **Metrics Calculation** (Stage 2)
   
   - Calculates **ARI** (Adjusted Rand Index) for each tool:
     - Measures similarity between predicted and true clone assignments
     - Range: -1 to 1, where 1 indicates perfect agreement
     - Accounts for chance agreement

4. **Visualization** (Stage 3)
   
   - **Confusion matrices**: Shows predicted vs. true clone assignments for each tool in a grid layout
   - **Bar plots**: ARI scores for each tool, with values displayed on top of bars

**Key Features**

- **Standardized input format**: All tools output predictions to TSV with ``barcode`` and ``prediction`` columns
- **Cell overlap tracking**: Ensures valid comparisons by using common cell set
- **Flexible clustering**: Tools that don't output clone labels can use hierarchical clustering with specified ``k``
- **Group-aware processing**: InferCNV processes multiple groups separately and combines results
- **Automated visualization**: Generates publication-quality figures for performance comparison
