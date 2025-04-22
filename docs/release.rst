..
   History
   =======



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
