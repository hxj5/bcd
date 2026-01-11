# cna_calling.py - CNA calling pipeline.


import anndata as ad
import gc
import numpy as np
import os
import pandas as pd
import shutil
import stat
import subprocess
import sys
import time
from logging import info, error
from logging import warning as warn


############################################
#------------------ main ------------------#
############################################

def cna_calling_main(
    sample_id,
    bam_fn,
    barcode_fn,
    cell_anno_fn,
    ref_cell_types,
    spot_pos_fn,
    out_dir,
    ncores,
    genome,
    umi_len,
    strand,
    star_genome_dir,
    infercnv_feature_fn,
    numbat_pp_script_fn,
    xcltk_feature_fn,
    eagle_fn,
    gmap_fn,
    snp_vcf_fn,
    panel_dir,
    fig_dpi = 300,
    script_dir = None,
    verbose = True
):
    os.makedirs(out_dir, exist_ok = True)

    if script_dir is None:
        script_dir = os.path.join(out_dir, "scripts")
    os.makedirs(script_dir, exist_ok = True)


    # run STARsolo
    if verbose:
        info("run STARsolo ...")

    rdr_dir = os.path.join(out_dir, "rdr_starsolo")
    os.makedirs(rdr_dir, exist_ok = True)

    star_cb_whitelist_fn = os.path.join(
        rdr_dir,
        "%s.starsolo.barcodes.pseudo_whitelist.tsv" % sample_id
    )
    with open(barcode_fn, "r") as fp:
        barcodes = [s.strip().split('-')[0] + '\n' for s in fp.readlines()]
    with open(star_cb_whitelist_fn, "w") as fp:
        fp.write("".join(barcodes))

    run_starsolo(
        bam_fn = bam_fn,
        genome_dir = star_genome_dir,
        cb_whitelist_fn = star_cb_whitelist_fn,
        out_dir = rdr_dir,
        ncores = ncores,
        read_file_type = 'SAM SE',
        read_file_command = 'samtools view -F 0x100',
        solo_strand = strand,
        solo_type = 'CB_UMI_Simple',
        solo_barcode_tags = 'CR UR',
        solo_qual_tags = 'CY UY',
        solo_umi_len = umi_len,
        out_sam_type = 'BAM Unsorted',
        script_fn = os.path.join(script_dir, "star.align.sh")
    )


    # run CopyKAT
    if verbose:
        info("run CopyKAT ...")

    copykat_dir = os.path.join(out_dir, "copykat")
    os.makedirs(copykat_dir, exist_ok = True)

    run_copykat(
        sample_id = sample_id,
        matrix_dir = rdr_dir,
        out_dir = copykat_dir,
        cell_anno_fn = cell_anno_fn,
        ref_cell_types = ref_cell_types,
        ncores = ncores,
        script_fn = os.path.join(script_dir, "copykat.call.R")
    )


    # Settings of inferCNV
    if verbose:
        info("run inferCNV ...")

    infercnv_dir = os.path.join(out_dir, "infercnv")
    os.makedirs(infercnv_dir, exist_ok = True)

    run_infercnv(
        sample_id = sample_id,
        matrix_dir = rdr_dir,
        out_dir = infercnv_dir,
        cell_anno_fn = cell_anno_fn,
        ref_cell_types = ref_cell_types,
        gene_anno_fn = infercnv_feature_fn,
        ncores = ncores,
        cutoff = 0.1,
        png_dpi = fig_dpi,
        script_fn = os.path.join(script_dir, "infercnv.call.R")
    )


    # Settings of Numbat
    if verbose:
        info("run Numbat ...")

    numbat_dir = os.path.join(out_dir, "numbat")
    os.makedirs(numbat_dir, exist_ok = True)

    run_numbat(
        sample_id = sample_id,
        bam_fn = bam_fn,
        barcode_fn = barcode_fn,
        matrix_dir = rdr_dir,
        out_dir = numbat_dir,
        cell_anno_fn = cell_anno_fn,
        ref_cell_types = ref_cell_types,
        ncores = ncores,
        pp_script_fn = numbat_pp_script_fn,
        eagle_fn = eagle_fn,
        gmap_fn = gmap_fn,
        snp_vcf_fn = snp_vcf_fn,
        panel_dir = panel_dir,
        platform = "10x",
        genome = genome,
        gamma = 20,
        script_dir = script_dir,
    )


    # Settings of XClone
    if verbose:
        info("run XClone ...")

    xclone_dir = os.path.join(out_dir, "xclone")
    os.makedirs(xclone_dir, exist_ok = True)

    xclone_cell_anno_fn = os.path.join(xclone_dir, "xclone.cell_anno.csv")
    df = pd.read_csv(cell_anno_fn, sep = '\t', header = None)
    df.columns = ['cell', 'cell_type']
    df.to_csv(xclone_cell_anno_fn, index = False, header = True)

    run_xclone(
        sample_id = sample_id,
        bam_fn = bam_fn,
        barcode_fn = barcode_fn,
        snp_vcf_fn = snp_vcf_fn,
        feature_fn = xcltk_feature_fn,
        cell_anno_fn = xclone_cell_anno_fn,
        ref_cell_types = ref_cell_types,
        out_dir = xclone_dir,
        ncores = ncores,
        eagle_fn = eagle_fn,
        gmap_fn = gmap_fn,
        panel_dir = panel_dir,
        genome = genome,
        cell_anno_key = "cell_type",
        plot_cell_anno_key = "cell_type",
        barcode_key = "cell",
        spot_pos_fn = spot_pos_fn,
        set_spatial = False if spot_pos_fn is None else True,
        xclone_plot = True,
        script_dir = script_dir
    )

    if verbose:
        info("All Done!")

    return(0)



################################################
#------------------ CalicoST ------------------#
################################################



###############################################
#------------------ CopyKAT ------------------#
###############################################

def run_copykat(
    sample_id,
    matrix_dir,
    out_dir,
    cell_anno_fn,
    ref_cell_types,
    ncores,
    script_fn = None
):
    s  = '''# run CopyKAT.\n'''
    s += '''\n'''
    s += '''library(copykat)\n'''
    s += '''library(Seurat)\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''run_copykat <- function(\n'''
    s += '''    sample_id,\n'''
    s += '''    matrix_dir, out_dir,\n'''
    s += '''    cell_anno_fn, ref_cell_types,\n'''
    s += '''    ncores)\n'''
    s += '''{\n'''
    s += '''    # check args\n'''
    s += '''    if (! dir.exists(out_dir)) {\n'''
    s += '''        dir.create(out_dir, recursive = T)\n'''
    s += '''    }\n'''
    s += '''    setwd(out_dir)\n'''
    s += '''    set.seed(123)\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''    # load matrix with Seurat\n'''
    s += '''    raw.data <- Read10X(data.dir = matrix_dir)\n'''
    s += '''    raw.data <- CreateSeuratObject(\n'''
    s += '''        counts = raw.data,\n'''
    s += '''        project = sample_id,\n'''
    s += '''        min.cells = 0,\n'''
    s += '''        min.features = 0\n'''
    s += '''    )\n'''
    s += '''    exp_raw_mtx <- tryCatch(\n'''
    s += '''        as.matrix(raw.data@assays$RNA@counts),\n'''
    s += '''        # work with Seurat v5.\n'''
    s += '''        # ref: https://github.com/satijalab/seurat/issues/9176\n'''
    s += '''        error = function(e) {\n'''
    s += '''            as.matrix(LayerData(\n'''
    s += '''                object = raw.data,\n'''
    s += '''                layer = "counts",\n'''
    s += '''                assay = "RNA")\n'''
    s += '''            )\n'''
    s += '''        }\n'''
    s += '''    )\n'''
    s += '''    saveRDS(exp_raw_mtx, file = paste0(out_dir, '/', sample_id, '.gex.rds'))\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''    # run copykat\n'''
    s += '''    ref_cells <- ""\n'''
    s += '''    if (! is.null(ref_cell_types)) {\n'''
    s += '''        cell_anno <- read.table(cell_anno_fn, sep = "\t", header = F, stringsAsFactors = F)\n'''
    s += '''        colnames(cell_anno) <- c("cell", "cell_type")\n'''
    s += '''        ref_cells <- cell_anno$cell[cell_anno$cell_type %in% ref_cell_types]\n'''
    s += '''    }\n'''
    s += '''    copykat_obj <- copykat(\n'''
    s += '''        rawmat = exp_raw_mtx,\n'''
    s += '''        id.type = "S",\n'''
    s += '''        ngene.chr = 5,\n'''
    s += '''        win.size = 25,\n'''
    s += '''        KS.cut = 0.1,\n'''
    s += '''        sam.name = sample_id,\n'''
    s += '''        distance = "euclidean",\n'''
    s += '''        norm.cell.names = ref_cells,\n'''
    s += '''        n.cores = ncores\n'''
    s += '''    )\n'''
    s += '''    saveRDS(copykat_obj, paste0(out_dir, '/', sample_id, '.copykat.obj.rds'))\n'''
    s += '''}\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''run_copykat(\n'''
    s += '''    sample_id = "%s",\n''' % sample_id
    s += '''    matrix_dir = "%s",\n''' % matrix_dir
    s += '''    out_dir = "%s",\n''' % out_dir
    s += '''    cell_anno_fn = "%s",\n''' % cell_anno_fn

    if ref_cell_types is None:
        s += '''    ref_cell_types = NULL,\n'''
    elif isinstance(ref_cell_types, str):
        s += '''    ref_cell_types = c("%s"),\n''' % ref_cell_types
    else:
        s += '''    ref_cell_types = %s,\n''' % \
            iter2c(ref_cell_types, dtype = 'str')

    s += '''    ncores = %d\n''' % ncores
    s += ''')\n'''
    s += '''\n'''
    s += '''print("[CopyKAT] All Done!")\n'''
    s += '''\n'''

    os.makedirs(out_dir, exist_ok = True)

    if script_fn is None:
        script_dir = os.path.join(out_dir, "scripts")
        os.makedirs(script_dir, exist_ok = True)
        script_fn = os.path.join(script_dir, "copykat.call.R")

    with open(script_fn, "w") as fp:
        fp.write(s)
    set_file_exe(script_fn)

    cmd = "/usr/bin/time -v Rscript %s" % script_fn
    ret, outs, errs = run_cmd(cmd)
    return(ret)



################################################
#------------------ InferCNV ------------------#
################################################

def run_infercnv(
    sample_id,
    matrix_dir,
    out_dir,
    cell_anno_fn,
    ref_cell_types,
    gene_anno_fn,
    ncores,
    cutoff = 0.1,
    png_dpi = 300,
    script_fn = None
):
    s  = '''# run inferCNV.\n'''
    s += '''\n'''
    s += '''library(infercnv)\n'''
    s += '''library(Seurat)\n'''
    s += '''\n'''
    s += '''run_infercnv <- function(\n'''
    s += '''    sample_id,\n'''
    s += '''    matrix_dir, out_dir,\n'''
    s += '''    cell_anno_fn, ref_cell_types, gene_anno_fn,\n'''
    s += '''    ncores,\n'''
    s += '''    cutoff, png_dpi)\n'''
    s += '''{\n'''
    s += '''    # check args\n'''
    s += '''    if (! dir.exists(out_dir)) {\n'''
    s += '''        dir.create(out_dir, recursive = T)\n'''
    s += '''    }\n'''
    s += '''    setwd(out_dir)\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''    # run infercnv\n'''
    s += '''    gex_mtx <- Seurat::Read10X(data.dir = matrix_dir)    # gene x cell\n'''
    s += '''    infercnv_obj <- CreateInfercnvObject(\n'''
    s += '''        raw_counts_matrix = gex_mtx,\n'''
    s += '''        annotations_file = cell_anno_fn,\n'''
    s += '''        delim = '\\t',\n'''
    s += '''        gene_order_file = gene_anno_fn,\n'''
    s += '''        ref_group_names = ref_cell_types\n'''
    s += '''    )\n'''
    s += '''    infercnv_obj <- infercnv::run(\n'''
    s += '''        infercnv_obj,\n'''
    s += '''        cutoff = cutoff,\n'''
    s += '''        out_dir = out_dir,\n'''
    s += '''        cluster_by_groups = T,\n'''
    s += '''        denoise = T,\n'''
    s += '''        HMM = T,\n'''
    s += '''        num_threads = ncores,\n'''
    s += '''        png_res = png_dpi\n'''
    s += '''    )\n'''
    s += '''}\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''run_infercnv(\n'''
    s += '''    sample_id = "%s",\n''' % sample_id
    s += '''    matrix_dir = "%s",\n''' % matrix_dir
    s += '''    out_dir = "%s",\n''' % out_dir
    s += '''    cell_anno_fn = "%s",\n''' % cell_anno_fn

    if ref_cell_types is None:
        s += '''    ref_cell_types = NULL,\n'''
    elif isinstance(ref_cell_types, str):
        s += '''    ref_cell_types = c("%s"),\n''' % ref_cell_types
    else:
        s += '''    ref_cell_types = %s,\n''' % \
            iter2c(ref_cell_types, dtype = 'str')

    s += '''    gene_anno_fn = "%s",\n''' % gene_anno_fn
    s += '''    ncores = %d,\n''' % ncores
    s += '''    cutoff = %f,\n''' % cutoff
    s += '''    png_dpi = %d\n''' % png_dpi
    s += ''')\n'''
    s += '''\n'''
    s += '''print("[inferCNV] All Done!")\n'''
    s += '''\n'''

    os.makedirs(out_dir, exist_ok = True)

    if script_fn is None:
        script_dir = os.path.join(out_dir, "scripts")
        os.makedirs(script_dir, exist_ok = True)
        script_fn = os.path.join(script_dir, "infercnv.call.R")

    with open(script_fn, "w") as fp:
        fp.write(s)
    set_file_exe(script_fn)

    cmd = "/usr/bin/time -v Rscript %s" % script_fn
    ret, outs, errs = run_cmd(cmd)
    return(ret)



##############################################
#------------------ Numbat ------------------#
##############################################

def run_numbat(
    sample_id,
    bam_fn,
    barcode_fn,
    matrix_dir,
    out_dir,
    cell_anno_fn,
    ref_cell_types,
    ncores,
    pp_script_fn,     # the pileup_and_phasing.R
    eagle_fn,
    gmap_fn,
    snp_vcf_fn,
    panel_dir,
    platform,
    genome,
    gamma = 20,
    script_dir = None
):
    os.makedirs(out_dir, exist_ok = True)

    pp_dir = os.path.join(out_dir, "allele")
    os.makedirs(pp_dir, exist_ok = True)

    if script_dir is None:
        script_dir = os.path.join(out_dir, "scripts")
    os.makedirs(script_dir, exist_ok = True)

    pp_script_fn_new = copy_file(pp_script_fn, script_dir)

    run_numbat_pp(
        sample_id = sample_id,
        bam_fn = bam_fn,
        barcode_fn = barcode_fn,
        out_dir = pp_dir,
        ncores = ncores,
        pp_script_fn = pp_script_fn_new,
        eagle_fn = eagle_fn,
        gmap_fn = gmap_fn,
        snp_vcf_fn = snp_vcf_fn,
        panel_dir = panel_dir,
        platform = platform,
        script_fn = os.path.join(script_dir, "numbat.pp.sh")
    )

    ret = run_numbat_call(
        sample_id = sample_id,
        allele_fn = os.path.join(pp_dir, "%s_allele_counts.tsv.gz" % sample_id),
        matrix_dir = matrix_dir,
        out_dir = os.path.join(out_dir, "cna"),
        cell_anno_fn = cell_anno_fn,
        ref_cell_types = ref_cell_types,
        ncores = ncores,
        out_prefix = "%s.numbat" % sample_id,
        genome = genome,
        gamma = gamma,
        script_fn = os.path.join(script_dir, "numbat.call.R")
    )
    return(ret)



def run_numbat_call(
    sample_id,
    allele_fn,
    matrix_dir,
    out_dir,
    cell_anno_fn,
    ref_cell_types,
    ncores,
    out_prefix,
    genome,
    gamma = 20,
    script_fn = None
):
    s  = '''# run Numbat CNA calling.\n'''
    s += '''\n'''
    s += '''library(numbat)\n'''
    s += '''library(Seurat)\n'''
    s += '''\n'''
    s += '''run_numbat_call <- function(\n'''
    s += '''    allele_fn, count_mtx_dir, out_dir,\n'''
    s += '''    cell_anno_fn, ref_cell_types,\n'''
    s += '''    ncores, out_prefix,\n'''
    s += '''    genome, gamma)\n'''
    s += '''{\n'''
    s += '''    # check args\n'''
    s += '''    if (! dir.exists(out_dir)) {\n'''
    s += '''        dir.create(out_dir, recursive = T)\n'''
    s += '''    }\n'''
    s += '''    setwd(out_dir)\n'''
    s += '''    if (! genome %in% c("hg19", "hg38"))\n'''
    s += '''        stop(sprintf("invalid genome version '%s'.", genome))\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''    # filter ref cells from allele dataframe\n'''
    s += '''    baf <- read.delim(allele_fn, header = T, stringsAsFactors = F)\n'''
    s += '''    cell_anno <- read.delim(cell_anno_fn, header = F, stringsAsFactors = F)\n'''
    s += '''    colnames(cell_anno) <- c("cell", "group")\n'''
    s += '''    ref_cells <- cell_anno[cell_anno$group %in% ref_cell_types, ]\n'''
    s += '''    baf_noref <- baf[! (baf$cell %in% ref_cells$cell), ]\n'''
    s += '''    saveRDS(baf_noref, sprintf('''
    s += '''        "%s/%s.ref_filtered.allele.df.rds", out_dir, out_prefix))\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''    # filter ref cells from count matrix\n'''
    s += '''    rdr <- Seurat::Read10X(data.dir = count_mtx_dir)   # gene x cell\n'''
    s += '''    rdr_noref <- rdr[, ! (colnames(rdr) %in% ref_cells$cell)]\n'''
    s += '''    saveRDS(rdr_noref, sprintf('''
    s += '''        "%s/%s.ref_filtered.count.mtx.rds", out_dir, out_prefix))\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''    # calculate average expression of ref cells\n'''
    s += '''    ref_mean_expr <- numbat::aggregate_counts(rdr, ref_cells)\n'''
    s += '''    saveRDS(ref_mean_expr, sprintf('''
    s += '''        "%s/%s.ref.gene_by_celltype.mtx.rds", out_dir, out_prefix))\n'''
    s += '''\n'''
    s += '''    obj <- run_numbat(\n'''
    s += '''        rdr_noref,         # gene x cell integer UMI count matrix \n'''
    s += '''        ref_mean_expr,     # reference expression profile, a gene x cell type normalized expression level matrix\n'''
    s += '''        baf_noref,         # allele dataframe generated by pileup_and_phase script\n'''
    s += '''        genome = genome,\n'''
    s += '''        t = 1e-5,\n'''
    s += '''        gamma = gamma,\n'''
    s += '''        ncores = ncores,\n'''
    s += '''        plot = TRUE,\n'''
    s += '''        out_dir = out_dir\n'''
    s += '''    )\n'''
    s += '''    saveRDS(obj, sprintf("%s/%s.out.object.rds", out_dir, out_prefix))\n'''
    s += '''}\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''run_numbat_call(\n'''
    s += '''    allele_fn = "%s",\n''' % allele_fn
    s += '''    count_mtx_dir = "%s",\n''' % matrix_dir
    s += '''    out_dir = "%s",\n''' % out_dir
    s += '''    cell_anno_fn = "%s",\n''' % cell_anno_fn

    if ref_cell_types is None:
        s += '''    ref_cell_types = NULL,\n'''
    elif isinstance(ref_cell_types, str):
        s += '''    ref_cell_types = c("%s"),\n''' % ref_cell_types
    else:
        s += '''    ref_cell_types = %s,\n''' % \
            iter2c(ref_cell_types, dtype = 'str')

    s += '''    ncores = %d,\n''' % ncores
    s += '''    out_prefix = "%s",\n''' % out_prefix
    s += '''    genome = "%s",\n''' % genome
    s += '''    gamma = %f\n''' % gamma
    s += ''')\n'''
    s += '''\n'''
    s += '''print("[Numbat] All Done!")\n'''
    s += '''\n'''

    os.makedirs(out_dir, exist_ok = True)

    if script_fn is None:
        script_fn = os.path.join(out_dir, "numbat.call.R")

    with open(script_fn, "w") as fp:
        fp.write(s)
    set_file_exe(script_fn)

    cmd = "/usr/bin/time -v Rscript %s" % script_fn
    ret, outs, errs = run_cmd(cmd)
    return(ret)



def run_numbat_pp(
    sample_id,
    bam_fn,
    barcode_fn,
    out_dir,
    ncores,
    pp_script_fn,
    eagle_fn,
    gmap_fn,
    snp_vcf_fn,
    panel_dir,
    platform,
    script_fn = None
):
    s  = '''#!/bin/bash\n'''
    s += '''# run Numbat preprocessing.\n'''
    s += '''\n'''
    s += '''Rscript  %s   \\\n''' % pp_script_fn
    s += '''    --label  %s    \\\n''' % sample_id
    s += '''    --samples  %s    \\\n''' % sample_id
    s += '''    --bams  %s    \\\n''' % bam_fn
    s += '''    --barcodes  %s    \\\n''' % barcode_fn
    s += '''    --outdir  %s    \\\n''' % out_dir
    s += '''    --ncores  %s    \\\n''' % ncores
    s += '''    --eagle  %s    \\\n''' % eagle_fn
    s += '''    --gmap  %s    \\\n''' % gmap_fn
    s += '''    --snpvcf  %s    \\\n''' % snp_vcf_fn
    s += '''    --paneldir  %s    \\\n''' % panel_dir

    if platform == "10x":
        s += "\n"
    elif platform == "smartseq":
        s += "--smartseq    \n"
    elif platform == "bulk":
        s += "--bulk    \n"
    else:
        raise ValueError("unknown platform '%s'!" % platform)

    s += '''\n'''

    os.makedirs(out_dir, exist_ok = True)

    if script_fn is None:
        script_fn = os.path.join(out_dir, "numbat.pp.sh")

    with open(script_fn, "w") as fp:
        fp.write(s)
    set_file_exe(script_fn)

    cmd = "/usr/bin/time -v %s" % script_fn
    ret, outs, errs = run_cmd(cmd)
    return(ret)



##############################################
#------------------ XClone ------------------#
##############################################

def run_xclone(
    sample_id,
    bam_fn,
    barcode_fn,
    snp_vcf_fn,
    feature_fn,
    cell_anno_fn,
    ref_cell_types,
    spot_pos_fn,
    out_dir,
    ncores,
    eagle_fn,
    gmap_fn,
    panel_dir,
    genome,
    cell_anno_key,
    plot_cell_anno_key,
    barcode_key,
    set_spatial = True,
    xclone_plot = True,
    script_dir = None
):
    os.makedirs(out_dir, exist_ok = True)

    baf_dir = os.path.join(out_dir, "xcltk_baf")
    os.makedirs(baf_dir, exist_ok = True)

    if script_dir is None:
        script_dir = os.path.join(out_dir, "scripts")
    os.makedirs(script_dir, exist_ok = True)

    if isinstance(ref_cell_types, str):
        ref_cell_types = [ref_cell_types]
    ref_cell_fn = os.path.join(baf_dir, "ref_cells.tsv")
    df = pd.read_csv(cell_anno_fn)
    df = df.loc[df[cell_anno_key].isin(ref_cell_types)].copy()
    df[[barcode_key]].to_csv(ref_cell_fn, index = False, header = False)


    ret = run_xcltk_baf(
        sample_id = sample_id,
        bam_fn = bam_fn,
        barcode_fn = barcode_fn,
        snp_vcf_fn = snp_vcf_fn,
        feature_fn = feature_fn,
        out_dir = baf_dir,
        ncores = ncores,
        ref_cell_fn = ref_cell_fn,
        eagle_fn = eagle_fn,
        gmap_fn = gmap_fn,
        panel_dir = panel_dir,
        script_fn = os.path.join(script_dir, "xcltk.baf.py")
    )

    rdr_dir = os.path.join(out_dir, "xcltk_rdr")
    os.makedirs(rdr_dir, exist_ok = True)

    ret = run_xcltk_rdr(
        bam_fn = bam_fn,
        barcode_fn = barcode_fn,
        feature_fn = feature_fn,
        out_dir = rdr_dir,
        ncores = ncores,
        script_fn = os.path.join(script_dir, "xcltk.rdr.py")
    )

    cna_dir = os.path.join(out_dir, "xclone_cna")
    os.makedirs(cna_dir, exist_ok = True)

    ret = run_xclone_call(
        sample_id = sample_id,
        baf_dir = os.path.join(baf_dir, "3_baf_fc"),
        rdr_dir = rdr_dir,
        cell_anno_fn = cell_anno_fn,
        ref_cell_types = ref_cell_types,
        out_dir = cna_dir,
        ncores = ncores,
        genome = genome,
        cell_anno_key = cell_anno_key,
        plot_cell_anno_key = plot_cell_anno_key,
        barcode_key = barcode_key,
        spot_pos_fn = spot_pos_fn,
        set_spatial = set_spatial,
        xclone_plot = xclone_plot,
        script_fn = os.path.join(script_dir, "xclone.call.py")
    )
    return(ret)



def run_xclone_call(
    sample_id,
    baf_dir,
    rdr_dir,
    cell_anno_fn,
    ref_cell_types,
    spot_pos_fn,
    out_dir,
    ncores,
    genome,
    cell_anno_key,
    plot_cell_anno_key,
    barcode_key,
    set_spatial = True,
    xclone_plot = True,
    script_fn = None
):
    s  = '''# run XClone CNA calling.\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''import anndata as ad\n'''
    s += '''import gc\n'''
    s += '''import os\n'''
    s += '''import scanpy as sc\n'''
    s += '''import sys\n'''
    s += '''import xclone\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''def run_xclone_call(\n'''
    s += '''    sample_id,\n'''
    s += '''    baf_dir, rdr_dir,\n'''
    s += '''    cell_anno_fn, ref_cell_types,\n'''
    s += '''    out_dir,\n'''
    s += '''    ncores,\n'''
    s += '''    genome,\n'''
    s += '''    cell_anno_key, plot_cell_anno_key, barcode_key,\n'''
    s += '''    spot_pos_fn, set_spatial,\n'''
    s += '''    xclone_plot\n'''
    s += '''):\n'''
    s += '''    # pre-check\n'''
    s += '''    xclone.pp.efficiency_preview()\n'''
    s += '''    xp_config = xclone.PreprocessingConfig(\n'''
    s += '''        dataset_name = sample_id,\n'''
    s += '''        module = "pre_check",\n'''
    s += '''        rdr_data_dir = rdr_dir,\n'''
    s += '''        baf_data_dir = baf_dir\n'''
    s += '''    )\n'''
    s += '''    xp_config.display()\n'''
    s += '''    xclone.pp.load_Xdata(module = "pre_check",  config_file = xp_config)\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''    # load data\n'''
    s += '''    ## rdr\n'''
    s += '''    xp_config = xclone.PreprocessingConfig(\n'''
    s += '''        dataset_name = sample_id,\n'''
    s += '''        module = "RDR",\n'''
    s += '''        barcodes_key = barcode_key,\n'''
    s += '''        set_spatial = set_spatial,\n'''
    s += '''        spot_position_file = spot_pos_fn,\n'''
    s += '''        rdr_data_dir = rdr_dir\n'''
    s += '''    )\n'''
    s += '''    xp_config.genome_mode = "%s_genes" % genome\n'''
    s += '''    xp_config.cell_anno_file = cell_anno_fn\n'''
    s += '''    xp_config.cell_anno_key = cell_anno_key\n'''
    s += '''    xp_config.display()\n'''
    s += '''    RDR_adata = xclone.pp.load_Xdata(module = "RDR", config_file = xp_config)\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''    ## baf\n'''
    s += '''    xp_config = xclone.PreprocessingConfig(\n'''
    s += '''        dataset_name = sample_id,\n'''
    s += '''        module = "BAF",\n'''
    s += '''        barcodes_key = barcode_key,\n'''
    s += '''        set_spatial = set_spatial,\n'''
    s += '''        spot_position_file = spot_pos_fn,\n'''
    s += '''        baf_data_dir = baf_dir\n'''
    s += '''    )\n'''
    s += '''    xp_config.genome_mode = "%s_genes" % genome\n'''
    s += '''    xp_config.cell_anno_file = cell_anno_fn\n'''
    s += '''    xp_config.cell_anno_key = cell_anno_key\n'''
    s += '''    xp_config.display()\n'''
    s += '''    BAF_adata = xclone.pp.load_Xdata(module = "BAF", config_file = xp_config)\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''    # run RDR module\n'''
    s += '''    xconfig = xclone.XCloneConfig(\n'''
    s += '''        dataset_name = sample_id,\n'''
    s += '''        module = "RDR",\n'''
    s += '''        set_spatial = set_spatial\n'''
    s += '''    )\n'''
    s += '''    xconfig.set_figure_params(xclone = True, fontsize = 18)\n'''
    s += '''    xconfig.outdir = out_dir\n'''
    s += '''    xconfig.cell_anno_key = cell_anno_key\n'''
    s += '''    xconfig.ref_celltype = ref_cell_types\n'''
    s += '''    xconfig.top_n_marker = 15 ##\n'''
    s += '''    xconfig.filter_ref_ave = 1.8 ##\n'''
    s += '''    xconfig.min_gene_keep_num = 1000 ##\n'''
    s += '''    xconfig.marker_group_anno_key = cell_anno_key\n'''
    s += '''    xconfig.xclone_plot = xclone_plot\n'''
    s += '''    xconfig.plot_cell_anno_key = plot_cell_anno_key\n'''
    s += '''    xconfig.fit_GLM_libratio = False\n'''
    s += '''    xconfig.exclude_XY = False\n'''
    s += '''    xconfig.remove_guide_XY = True\n'''
    s += '''    xconfig.display()\n'''
    s += '''    RDR_Xdata = xclone.model.run_RDR(RDR_adata, config_file = xconfig)\n'''
    s += '''\n'''
    s += '''    # to save memory as following BAF module will use multiprocessing.\n'''
    s += '''    del RDR_Xdata\n'''
    s += '''    gc.collect()\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''    # run BAF module\n'''
    s += '''    xconfig = xclone.XCloneConfig(\n'''
    s += '''        dataset_name = sample_id,\n'''
    s += '''        module = "BAF",\n'''
    s += '''        set_spatial = set_spatial\n'''
    s += '''    )\n'''
    s += '''    xconfig.update_info_from_rdr = False      # update using only BAF information.\n'''
    s += '''    xconfig.set_figure_params(xclone = True, fontsize = 18)\n'''
    s += '''    xconfig.outdir = out_dir\n'''
    s += '''    xconfig.cell_anno_key = cell_anno_key\n'''
    s += '''    xconfig.ref_celltype = ref_cell_types\n'''
    s += '''    xconfig.xclone_plot = xclone_plot\n'''
    s += '''    xconfig.plot_cell_anno_key = plot_cell_anno_key\n'''
    s += '''    xconfig.bin_nproc = ncores\n'''
    s += '''    xconfig.HMM_nproc = ncores\n'''
    s += '''    xconfig.display()\n'''
    s += '''    BAF_merge_Xdata = xclone.model.run_BAF(BAF_adata, config_file = xconfig)\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''    # run Combine module\n'''
    s += '''    xconfig = xclone.XCloneConfig(\n'''
    s += '''        dataset_name = sample_id,\n'''
    s += '''        module = "Combine",\n'''
    s += '''        set_spatial = set_spatial\n'''
    s += '''    )\n'''
    s += '''    xconfig.set_figure_params(xclone = True, fontsize = 18)\n'''
    s += '''    xconfig.outdir = out_dir\n'''
    s += '''    xconfig.cell_anno_key = cell_anno_key\n'''
    s += '''    xconfig.ref_celltype = ref_cell_types\n'''
    s += '''    xconfig.xclone_plot = xclone_plot\n'''
    s += '''    xconfig.plot_cell_anno_key = plot_cell_anno_key\n'''
    s += '''    xconfig.merge_loss = False\n'''
    s += '''    xconfig.merge_loh = False\n'''
    s += '''    xconfig.BAF_denoise = True\n'''
    s += '''    xconfig.exclude_XY = True\n'''
    s += '''    xconfig.display()\n'''
    s += '''\n'''
    s += '''    fn = os.path.join(out_dir, "data/RDR_adata_KNN_HMM_post.h5ad")\n'''
    s += '''    RDR_Xdata = ad.read_h5ad(fn)\n'''
    s += '''\n'''
    s += '''    combine_Xdata = xclone.model.run_combine(\n'''
    s += '''        RDR_Xdata,\n'''
    s += '''        BAF_merge_Xdata,\n'''
    s += '''        verbose = True,\n'''
    s += '''        run_verbose = True,\n'''
    s += '''        config_file = xconfig\n'''
    s += '''    )\n'''
    s += '''\n'''
    s += '''\n'''
    s += '''run_xclone_call(\n'''
    s += '''    sample_id = "%s",\n''' % sample_id
    s += '''    baf_dir = "%s/",\n''' % baf_dir
    s += '''    rdr_dir = "%s/",\n''' % rdr_dir
    s += '''    cell_anno_fn = "%s",\n''' % cell_anno_fn

    if isinstance(ref_cell_types, str):
        s += '''    ref_cell_types = "%s",\n''' % ref_cell_types
    else:
        s += '''    ref_cell_types = %s,\n''' % str(ref_cell_types)

    s += '''    out_dir = "%s/",\n''' % out_dir
    s += '''    ncores = %d,\n''' % ncores
    s += '''    genome = "%s",\n''' % genome
    s += '''    cell_anno_key = "%s",\n''' % cell_anno_key
    s += '''    plot_cell_anno_key = "%s",\n''' % plot_cell_anno_key
    s += '''    barcode_key = "%s",\n''' % barcode_key
    s += '''    spot_pos_fn = %s,\n''' % str_or_none(spot_pos_fn)
    s += '''    set_spatial = %s,\n''' % str(set_spatial)
    s += '''    xclone_plot = %s\n''' % str(xclone_plot)
    s += ''')\n'''
    s += '''\n'''
    s += '''print("[XClone] All Done!")\n'''
    s += '''\n'''

    os.makedirs(out_dir, exist_ok = True)

    if script_fn is None:
        script_fn = os.path.join(out_dir, "xclone.call.py")

    with open(script_fn, "w") as fp:
        fp.write(s)
    set_file_exe(script_fn)

    cmd = "/usr/bin/time -v python %s" % script_fn
    ret, outs, errs = run_cmd(cmd)
    return(ret)



def run_xcltk_baf(
    sample_id,
    bam_fn,
    barcode_fn,
    snp_vcf_fn,
    feature_fn,
    out_dir,
    ncores,
    ref_cell_fn,
    eagle_fn,
    gmap_fn,
    panel_dir,
    min_count = 11,
    min_maf = 0.1,
    script_fn = None
):
    cell_tag = "CB"
    umi_tag = "UB"

    s  = '''# run xcltk baf.\n'''
    s += '''\n'''
    s += '''from xcltk.baf.pipeline import pipeline_wrapper\n'''
    s += '''\n'''
    s += '''pipeline_wrapper(\n'''
    s += '''    label = "%s",\n''' % sample_id
    s += '''    sam_fn = "%s",\n''' % bam_fn
    s += '''    barcode_fn = "%s",\n''' % barcode_fn
    s += '''    snp_vcf_fn = "%s",\n''' % snp_vcf_fn
    s += '''    region_fn = "%s",\n''' % feature_fn
    s += '''    out_dir = "%s",\n''' % out_dir
    s += '''    ncores = %d,\n''' % ncores
    s += '''    ref_cell_fn = %s,\n''' % str_or_none(ref_cell_fn)
    s += '''    gmap_fn = "%s",\n''' % gmap_fn
    s += '''    eagle_fn = "%s",\n''' % eagle_fn
    s += '''    panel_dir = "%s",\n''' % panel_dir
    s += '''    cell_tag = %s,\n''' % str_or_none(cell_tag)
    s += '''    umi_tag = %s,\n''' % str_or_none(umi_tag)
    s += '''    min_count = %d,\n''' % min_count
    s += '''    min_maf = %f\n''' % min_maf
    s += ''')\n'''
    s += '''\n'''

    os.makedirs(out_dir, exist_ok = True)

    if script_fn is None:
        script_fn = os.path.join(out_dir, "xcltk.baf.py")

    with open(script_fn, "w") as fp:
        fp.write(s)
    set_file_exe(script_fn)

    cmd = "/usr/bin/time -v python %s" % script_fn
    ret, outs, errs = run_cmd(cmd)
    return(ret)



def run_xcltk_rdr(
    bam_fn,
    barcode_fn,
    feature_fn,
    out_dir,
    ncores,
    script_fn = None
):
    s  = '''# run xcltk rdr.\n'''
    s += '''\n'''
    s += '''from xcltk.rdr.fc.main import fc_wrapper\n'''
    s += '''\n'''
    s += '''fc_wrapper(\n'''
    s += '''    sam_fn = "%s",\n''' % bam_fn
    s += '''    barcode_fn = "%s",\n''' % barcode_fn
    s += '''    region_fn = "%s",\n''' % feature_fn
    s += '''    out_dir = "%s",\n''' % out_dir
    s += '''    ncores = %d\n''' % ncores
    s += ''')\n'''
    s += '''\n'''

    os.makedirs(out_dir, exist_ok = True)

    if script_fn is None:
        script_fn = os.path.join(out_dir, "xcltk.rdr.py")

    with open(script_fn, "w") as fp:
        fp.write(s)
    set_file_exe(script_fn)

    cmd = "/usr/bin/time -v python %s" % script_fn
    ret, outs, errs = run_cmd(cmd)
    return(ret)



################################################
#------------------ STARsolo ------------------#
################################################

def run_starsolo(
    bam_fn,
    genome_dir,
    cb_whitelist_fn,
    out_dir,
    ncores,
    read_file_type = 'SAM SE',
    read_file_command = 'samtools view -F 0x100',
    solo_strand = 'Forward',
    solo_type = 'CB_UMI_Simple',
    solo_barcode_tags = 'CR UR',
    solo_qual_tags = 'CY UY',
    solo_umi_len = 12,
    out_sam_type = 'BAM Unsorted',
    script_fn = None
):
    """
    cb_whiltelist_fn: no "-1" suffix.
    """
    os.makedirs(out_dir, exist_ok = True)

    res_dir = os.path.join(out_dir, "starsolo_out")
    os.makedirs(res_dir, exist_ok = True)

    s  = '''# run STAR align.\n'''
    s += '''\n'''
    s += '''STAR  \\\n'''
    s += '''    --runThreadN  %d  \\\n''' % ncores
    s += '''    --outFileNamePrefix  %s/  \\\n''' % res_dir
    s += '''    --genomeDir  %s  \\\n''' % genome_dir
    s += '''    --soloStrand  %s   \\\n''' % solo_strand
    s += '''    --soloType  %s  \\\n''' % solo_type
    s += '''    --soloCBwhitelist  %s  \\\n''' % cb_whitelist_fn
    s += '''    --readFilesIn  %s  \\\n''' % bam_fn
    s += '''    --readFilesType  %s  \\\n''' % read_file_type
    s += '''    --readFilesCommand  %s  \\\n''' % read_file_command
    s += '''    --soloInputSAMattrBarcodeSeq  %s  \\\n''' % solo_barcode_tags
    s += '''    --soloInputSAMattrBarcodeQual  %s  \\\n''' % solo_qual_tags
    s += '''    --soloUMIlen  %d    \\\n''' % solo_umi_len
    s += '''    --outSAMtype  %s    \n''' % out_sam_type
    s += '''\n'''

    if script_fn is None:
        script_dir = os.path.join(out_dir, "scripts")
        os.makedirs(script_dir, exist_ok = True)
        script_fn = os.path.join(script_dir, "star.align.sh")

    with open(script_fn, "w") as fp:
        fp.write(s)
    set_file_exe(script_fn)

    cmd = "/usr/bin/time -v %s" % script_fn
    ret, outs, errs = run_cmd(cmd)


    # the "Seurat::Read10X()" expects "genes.tsv" or "barcodes.tsv.gz".
    # "Seurat::Read10X()" is used by inferCNV and Numbat for loading matrix data.
    mtx_dir = os.path.join(res_dir, "Solo.out/Gene/raw")

    with open(os.path.join(mtx_dir, "barcodes.tsv"), "r") as fp:
        barcodes = [s.strip() + '-1\n' for s in fp.readlines()]
    with open(os.path.join(out_dir, "barcodes.tsv"), "w") as fp:
        fp.write("".join(barcodes))

    fn = os.path.join(mtx_dir, "features.tsv")
    df = pd.read_csv(fn, sep = '\t', header = None)
    df = df.iloc[:, [0, 1]].copy()     # df.columns = ['gene_id', 'gene_name']
    out_fn = os.path.join(out_dir, 'genes.tsv')
    df.to_csv(out_fn, sep = '\t', index = False, header = False)

    copy_file(os.path.join(mtx_dir, 'matrix.mtx'), out_dir)

    return(ret)



def star_index(
    fasta_fn,
    gtf_fn,
    out_dir,
    read_len,
    ncores,
    script_fn = None
):
    s  = '''# run STAR index.\n'''
    s += '''\n'''
    s += '''STAR  \\\n'''
    s += '''    --runThreadN  %d  \\\n''' % ncores
    s += '''    --runMode  genomeGenerate  \\\n'''
    s += '''    --genomeDir  %s  \\\n''' % out_dir
    s += '''    --genomeFastaFiles  %s  \\\n''' % fasta_fn
    s += '''    --sjdbGTFfile  %s  \\\n''' % gtf_fn
    s += '''    --sjdbOverhang  %d    \n''' % (read_len - 1, )
    s += '''\n'''

    os.makedirs(out_dir, exist_ok = True)

    if script_fn is None:
        script_dir = os.path.join(out_dir, "scripts")
        os.makedirs(script_dir, exist_ok = True)
        script_fn = os.path.join(script_dir, "star.index.sh")

    with open(script_fn, "w") as fp:
        fp.write(s)
    set_file_exe(script_fn)

    cmd = "/usr/bin/time -v %s" % script_fn
    ret, outs, errs = run_cmd(cmd)
    return(ret)



##################################################
#------------------ utils.py2r ------------------#
##################################################

def iter2c(items, dtype = 'int'):
    s = "c("
    for i, x in enumerate(items):
        if i > 0:
            s += ", "
        if dtype == 'int':
            s += "%d" % x
        elif dtype == 'str':
            s += '''"%s"''' % x
        else:
            raise ValueError()
    s += ")"
    return(s)



##################################################
#------------------ utils.py2r ------------------#
##################################################

def copy_file(src, dst):
    return shutil.copy2(src, dst)



def run_cmd(args, check_retcode = False):
    proc = subprocess.Popen(
        args = args,
        shell = True,
        executable = "/bin/bash",
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )
    outs, errs = proc.communicate()
    ret = proc.returncode
    if ret != 0:
        error(str(errs.decode()))
        if check_retcode:
            raise RuntimeError
    return(ret, outs, errs)



def set_file_exe(fn):
    st = os.stat(fn)
    os.chmod(fn, st.st_mode | stat.S_IXUSR)



def str_or_none(s):
    if s is None:
        return "None"
    else:
        return '"%s"' % str(s)
