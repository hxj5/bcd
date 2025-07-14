bcd - Benchmarking of CNA Detection from Single-cell and Spatial Transcriptomics
================================================================================
The ``bcd`` (Benchmarking of CNA Detection) pipeline evaluates the performance
of tools in detecting copy number alterations (CNAs) from single-cell and 
spatial transcriptomics, using the ROC and PRC metrics, 
given the input ground truth of CNA profiles.



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


Quick Usage
~~~~~~~~~~~

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



Acknowledgement
---------------
Previously, we have an `R version of the pipeline <https://github.com/Rongtingting/CNV_calling_Benchmark/tree/main/scripts/evaluate>`_.



.. _conda: https://docs.conda.io/en/latest/
.. _infercnv: https://github.com/broadinstitute/infercnv
