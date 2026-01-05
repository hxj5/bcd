# utils.py



# ref: https://github.com/kharchenkolab/numbat/blob/main/R/vis.R#L44
def plot_psbulk_baf(
    df_snp,
    df_genome,
    ...
):
    """Visualize CNA BAF signal in pseudobulks, aggregating cells by clone.
    
    Inputs
    ------
    SNP data - A pandas.DataFrame with columns:
        chrom : str
            Chromosome.
        pos : int
            SNP position (1-based).
        hap : {0, 1}
            Haplotype index.
        AF : float
            Raw allele frequency (before phasing), i.e., AD/DP.
    Genome data - A pandas.DataFrame with columns:
        chrom : str
            Chromosome.
        arm : {'p', 'q'}
            Chromosome arm.
        length : int
            Length of the chromosome arm.
    """
    pass



def plot_psbulk_rdr(
    df_gene,
    df_genome,
    ...
):
    """Visualize CNA RDR signal in pseudobulks, aggregating cells by clone.
    
    Inputs
    ------
    Gene data - A pandas.DataFrame with columns:
        chrom : str
            Chromosome.
        start : int
            Gene start position (1-based, inclusive).
        end : int
            Gene end position (1-based, inclusive).
        gene : str
            Gene name.
        log2FC : float
            Log2 fold change compared to reference cells.
    Genome data - A pandas.DataFrame with columns:
        chrom : str
            Chromosome.
        arm : {'p', 'q'}
            Chromosome arm.
        length : int
            Length of the chromosome arm.
    """
    pass



def plot_psbulk(
    df_gene,
    df_snp,
    df_genome,
    ...
):
    pass
 