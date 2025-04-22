# gscale.py - convert between different genomic scales, e.g., region and gene.


from logging import info
from logging import warning as warn
from .grange import format_chrom



def get_overlap_genes(df, anno):
    """Get overlapping genes of regions.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Region-scale data.
        It should contain at least four columns:
        - "chrom";
        - "start": 1-based and inclusive;
        - "end": 1-based and inclusive;
        - "region": region ID.
    anno : pandas.DataFrame
        Gene annotations.
        It should contain at least four columns:
        - "chrom";
        - "start": 1-based and inclusive;
        - "end": 1-based and inclusive;
        - "gene": gene name.
        
    Returns
    -------
    dict
        Overlapping results.
    """
    df["chrom"] = df["chrom"].map(format_chrom)
    anno["chrom"] = anno["chrom"].map(format_chrom)
    
    df = df.drop_duplicates("region", ignore_index = True)
    anno = anno.drop_duplicates("gene", ignore_index = True)
    
    overlap = df.groupby("region").apply(
        lambda x: anno.loc[
            (anno["chrom"] == x["chrom"].values[0]) &
            (anno["start"] <= x["end"].values[0]) &
            (anno["end"] >= x["start"].values[0]),
            "gene"
        ]).reset_index()
    overlap = overlap[["region", "gene"]].sort_values(by = "region")
    
    stat = overlap.groupby("gene").size().reset_index(name = "n")
    dup = stat[stat["n"] > 1]
    uniq = stat[stat["n"] == 1].merge(overlap, on = "gene", how = "left")
    uniq = uniq[["region", "gene"]].sort_values(by = "region")
    
    res = dict(
        # overlap : pandas.DataFrame
        #   The overlapping results. It contains two columns:
        #   - "region": region ID.
        #   - "gene": name of genes overlapping the region.
        overlap = overlap,
        
        # n_region : int
        #   Number of unique regions.
        n_region = df.shape[0],
        
        # n_region_overlap : int
        #   Number of unique regions that have overlapping genes.
        n_region_overlap = len(overlap["region"].unique()),
        
        # n_gene : int
        #   Number of unique genes.
        n_gene = anno.shape[0],
        
        # n_gene_overlap : int
        #   Number of unique genes that have overlapping regions.
        n_gene_overlap = len(overlap["gene"].unique()),
        
        # n_gene_dup : int
        #   Number of genes overlapping more than 1 regions.
        n_gene_dup = dup.shape[0],
        
        # overlap_uniq : pandas.DataFrame
        #   The overlapping results. 
        #   Similar to `overlap`, but the genes overlapping more than 1 regions
        #   are removed.
        overlap_uniq = uniq
    )
    
    return(res)



def reg2gene(df, anno, no_dup = True, verbose = True):
    """Convert non-overlapping region-scale data into gene-scale.
    
    Each row/record in `df` will be copied (and updated accordingly) for
    each gene overlapping with the "region" within the record.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Non-overlapping region-scale data, i.e., the regions should not 
        overlap each other.
        It should contain at least four columns:
        - "chrom";
        - "start": 1-based and inclusive;
        - "end": 1-based and inclusive;
        - "region": region ID.
    anno : pandas.DataFrame
        Gene annotations.
        It should contain at least four columns:
        - "chrom";
        - "start": 1-based and inclusive;
        - "end": 1-based and inclusive;
        - "gene": gene name.
    no_dup : bool, default True
        Whether to remove genes overlapping more than 1 regions from final 
        result.
    verbose : bool, default True
        Whether to show detailed logging information.
        
    Returns
    -------
    dict
        The convertted results.
    """
    df["chrom"] = df["chrom"].map(format_chrom)
    anno["chrom"] = anno["chrom"].map(format_chrom)
    
    res = get_overlap_genes(df, anno)
    overlap = res["overlap"]
    if no_dup:
         overlap = res["overlap_uniq"]
    
    s = "\n"
    if verbose:
        s += "\tsettings: no_dup = %s;\n" % str(no_dup)
        s += "\t#regions: unique = %d;\n" % res["n_region"]
        s += "\t#regions: having overlap genes = %d;\n" % res["n_region_overlap"]
        s += "\t#genes: unique = %d;\n" % res["n_gene"]
        s += "\t#genes: having overlap regions = %d;\n" % res["n_gene_overlap"]
        s += "\t#genes: overlap more than 1 regions = %d;\n" % res["n_gene_dup"]
        s += "\t#genes: overlap exactly 1 regions = %d;\n" % res["overlap_uniq"].shape[0]
        s += "\t#records: input = %d;\n" % df.shape[0]
        
        
    df = df.merge(overlap, on = "region", how = "inner")
    
    if verbose:    
        s += "\t#records: output = %d;\n" % df.shape[0]
        s += "\t#regions: output unique = %d;\n" % len(df["region"].unique())
        s += "\t#genes: output unique = %d;\n" % len(df["gene"].unique())
        
        info("reg2gene results:\n%s" % s)
    
    
    res = dict(
        # df : pandas.DataFrame
        #   The converted data. It contains at least five columns:
        #   - "chrom";
        #   - "start": start pos of the region; 1-based and inclusive;
        #   - "end": end pos of the region; 1-based and inclusive;
        #   - "region": region ID.
        #   - "gene": name of gene overlapping the region.
        df = df,
        
        # overlap : pandas.DataFrame
        #   The overlapping results. It contains two columns:
        #   - "region": region ID.
        #   - "gene": name of genes overlapping the region. 
        overlap = overlap
    )

    return(res)
