..
   History
   =======
   

Release v0.5.0 (29/12/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
For each module, merge all files into one single file.


 
Release v0.4.0 (26/12/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Feature enhancement:

* make v0.3.0 a module named cna_profile;
* add another two modules: tumor_nontumor and subclonal_structure.


 
Release v0.3.0 (22/05/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Feature enhancement:

* add option ``overlap_how`` indicating how to subset the tool matrices given
  the overlap cells and genes.
  Default is "isec-cells": subset tool matrix by intersected cells only.
  Previously, subset tool matrix by intersected cells and genes (the 
  "isec-both" option value).

Others

* docs: update manual, adding options "numbat_mtx_how" and "overlap_how".
* restructure folder src_deprecated.



Release v0.2.0 (03/05/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Feature enhancement:

* add option ``numbat_mtx_how`` indicating how to process the extracted 
  Numbat matrix before overlap step.
  Default is "expand": expand the Numbat matrix to transcriptomics scale and 
  fill value 0.
  Previously, the raw Numbat matrix is used (the "raw" option value).

Fix bug:

* allow empty input truth file.
* extract: remove duplicate records in Numbat object file.
* overlap: skip all next steps if no any tools support some CNA type.
* utils.gscale.get_overlap_genes: select two columns to force returning 
  DataFrame instead of Series when there is only one region.

Others

* README: add dependency.



Release v0.1.0 (22/04/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This version restrcuture the framework, splitting the whole pipeline into
five steps, each wrapped in one module, including 
``extract, truth, overlap, metric, plot``.

* args: the input data of each tool should be wrapped into corresponding
  ``ToolArgs`` object.
* overlap: use intersected cells and genes of all tools and ground truth.
* metric: re-implement and update the R version of ``binaryROC``.

Note

* this version only supports inputs from inferCNV and Numbat.



Release v0.0.1 (16/04/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
* init converting R functions into python version.
