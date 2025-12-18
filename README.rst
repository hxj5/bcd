bcd - Benchmarking of CNA Detection from Single-cell and Spatial Transcriptomics
================================================================================
The ``bcd`` (Benchmarking of CNA Detection) pipeline evaluates the performance
of SOTA tools in CNA profile detection, tumor vs. non-tumor classification,
and subclonal structure identification.



News
----
Release notes are at `docs/release.rst <./docs/release.rst>`_.



Installation
------------

Dependency
~~~~~~~~~~
The pipeline depends on some tools and packages listed below:

* R (and Rscript), for extracting expression or probability matrix from 
  output of R packages, such as inferCNV and Numbat.


Install from this Github Repo (latest stable/dev version)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    pip install -U git+https://github.com/hxj5/bcd.git


Install locally
~~~~~~~~~~~~~~~

.. code-block:: shell

    git clone https://github.com/hxj5/bcd.git
    cd bcd
    pip install -U .


In either case, if you don't have write permission for your current Python
environment, we suggest creating a separate conda_ environment 
or add ``--user`` for your current one.



Manual
------
The full manual is at `docs/manual.rst <./docs/manual.rst>`_.



Acknowledgement
---------------
Previously, we have an `R version of the pipeline <https://github.com/Rongtingting/CNV_calling_Benchmark/tree/main/scripts/evaluate>`_.



.. _conda: https://docs.conda.io/en/latest/
.. _infercnv: https://github.com/broadinstitute/infercnv
