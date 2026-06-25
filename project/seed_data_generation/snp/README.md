
The list of phased SNPs (the `HCC3N.phased.vcf.gz` file) was generated
by the reference phasing pipeline (`scripts/pipeline.py`), 
mainly using cellsnp-lite for genotyping and Eagle2 for phasing that wrapped
in `xcltk.baf.genotype` module.

Specifically, it uses 4289 spots from the HCC3-N spatial transcriptomics
dataset for genotyping, with parameters `--minCOUNT 11` and `--minMAF 0.1`.

When running this pipeline with different parameters:

- mc10maf0.1: `--minCOUNT 10` and `--minMAF 0.1`; output 9665 phased SNPs.
- mc10maf0.2: `--minCOUNT 10` and `--minMAF 0.2`; output 8886 phased SNPs.
- mc11maf0.1: `--minCOUNT 11` and `--minMAF 0.1`; output 8866 phased SNPs.
- mc20maf0.1: `--minCOUNT 20` and `--minMAF 0.1`, output 5462 phased SNPs.

