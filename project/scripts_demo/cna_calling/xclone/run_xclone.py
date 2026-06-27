# run_xclone.py


import anndata as ad
import gc
import os
import pandas as pd
import shutil
import stat
import subprocess
import sys
import time
import xclone
from logging import info, error



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
    xclone_n_clones = None,
    set_spatial = False,
    xclone_plot = True,
    script_dir = None
):
    os.makedirs(out_dir, exist_ok = True)
    
    xclone_cell_anno_fn = os.path.join(out_dir, "xclone.cell_anno.csv")
    df = pd.read_csv(cell_anno_fn, sep = '\t', header = None)
    df.columns = ['cell', 'cell_type']
    df.to_csv(xclone_cell_anno_fn, index = False, header = True)
    cell_anno_fn = xclone_cell_anno_fn

    baf_dir = os.path.join(out_dir, "xcltk_baf/")
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

    rdr_dir = os.path.join(out_dir, "xcltk_rdr/")
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
        baf_dir = os.path.join(baf_dir, "3_baf_fc/"),
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
        xclone_n_clones = xclone_n_clones,
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
    xclone_n_clones = None,
    set_spatial = False,
    xclone_plot = True,
    script_fn = None
):
    os.makedirs(out_dir, exist_ok = True)

    try:
        # pre-check
        xclone.pp.efficiency_preview()
        xp_config = xclone.PreprocessingConfig(
            dataset_name = sample_id,
            module = "pre_check",
            rdr_data_dir = rdr_dir,
            baf_data_dir = baf_dir
        )
        xp_config.display()
        xclone.pp.load_Xdata(module = "pre_check", config_file = xp_config)

        # load data
        ## rdr
        xp_config = xclone.PreprocessingConfig(
            dataset_name = sample_id,
            module = "RDR",
            barcodes_key = barcode_key,
            set_spatial = set_spatial,
            spot_position_file = spot_pos_fn,
            rdr_data_dir = rdr_dir
        )
        xp_config.genome_mode = "%s_genes" % genome
        xp_config.cell_anno_file = cell_anno_fn
        xp_config.cell_anno_key = cell_anno_key
        xp_config.exclude_XY = True
        xp_config.display()
        RDR_adata = xclone.pp.load_Xdata(module = "RDR", config_file = xp_config)

        ## baf
        xp_config = xclone.PreprocessingConfig(
            dataset_name = sample_id,
            module = "BAF",
            barcodes_key = barcode_key,
            set_spatial = set_spatial,
            spot_position_file = spot_pos_fn,
            baf_data_dir = baf_dir
        )
        xp_config.genome_mode = "%s_genes" % genome
        xp_config.cell_anno_file = cell_anno_fn
        xp_config.cell_anno_key = cell_anno_key
        xp_config.exclude_XY = True
        xp_config.display()
        BAF_adata = xclone.pp.load_Xdata(module = "BAF", config_file = xp_config)

        # run RDR module
        xconfig = xclone.XCloneConfig(
            dataset_name = sample_id,
            module = "RDR",
            set_spatial = set_spatial
        )
        xconfig.set_figure_params(xclone = True, fontsize = 18)
        xconfig.outdir = out_dir
        xconfig.cell_anno_key = cell_anno_key
        xconfig.ref_celltype = ref_cell_types
        xconfig.marker_group_anno_key = cell_anno_key
        xconfig.xclone_plot = xclone_plot
        xconfig.plot_cell_anno_key = plot_cell_anno_key
        xconfig.n_clones = xclone_n_clones
        xconfig.exclude_XY = True
        xconfig.filter_ref_ave = 1.8

        xconfig.display()
        RDR_Xdata = xclone.model.run_RDR(RDR_adata, config_file = xconfig)

        # to save memory as following BAF module will use multiprocessing.
        del RDR_Xdata
        gc.collect()

        # run BAF module
        xconfig = xclone.XCloneConfig(
            dataset_name = sample_id,
            module = "BAF",
            set_spatial = set_spatial
        )
        xconfig.update_info_from_rdr = False      # update using only BAF information.
        xconfig.set_figure_params(xclone = True, fontsize = 18)
        xconfig.outdir = out_dir
        xconfig.cell_anno_key = cell_anno_key
        xconfig.ref_celltype = ref_cell_types
        xconfig.xclone_plot = xclone_plot
        xconfig.plot_cell_anno_key = plot_cell_anno_key
        xconfig.bin_nproc = ncores
        xconfig.HMM_nproc = ncores
        xconfig.exclude_XY = True
        xconfig.display()
        BAF_merge_Xdata = xclone.model.run_BAF(BAF_adata, config_file = xconfig)

        # run Combine module
        xconfig = xclone.XCloneConfig(
            dataset_name = sample_id,
            module = "Combine",
            set_spatial = set_spatial
        )
        xconfig.set_figure_params(xclone = True, fontsize = 18)
        xconfig.outdir = out_dir
        xconfig.cell_anno_key = cell_anno_key
        xconfig.ref_celltype = ref_cell_types
        xconfig.xclone_plot = xclone_plot
        xconfig.plot_cell_anno_key = plot_cell_anno_key
        xconfig.exclude_XY = True
        xconfig.display()

        fn = os.path.join(out_dir, "data/RDR_adata_KNN_HMM_post.h5ad")
        RDR_Xdata = ad.read_h5ad(fn)

        combine_Xdata = xclone.model.run_combine(
            RDR_Xdata,
            BAF_merge_Xdata,
            verbose = True,
            run_verbose = True,
            config_file = xconfig
        )

        info("[XClone] All Done!")
        return(0)

    except Exception as e:
        error("[XClone] Error occurred: %s" % str(e))
        return(1)



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



if __name__ == "__main__":
    run_xclone(
        sample_id = 'HCC3N_simu',
        bam_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/downsample_coverage/ds_50perc/gen_data/subset_bam/possorted.bam',
        barcode_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/downsample_coverage/ds_50perc/gen_data/rdr_starsolo/barcodes.tsv',
        snp_vcf_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/refapp/phasing/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz',
        feature_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/refapp/xcltk/xcltk.genes.tsv',
        cell_anno_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/base/data/matrix/spot_anno.tsv',
        ref_cell_types = 'normal',
        spot_pos_fn = None,
        out_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/downsample_coverage/ds_50perc/cna_calling/xclone',
        ncores = 10,
        eagle_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/refapp/eagle/eagle',
        gmap_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/refapp/eagle/genetic_map_hg38_withX.txt.gz',
        panel_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/refapp/phasing/1000G_hg38',
        genome = 'hg38',
        cell_anno_key = 'cell_type',
        plot_cell_anno_key = 'cell_type',
        barcode_key = 'cell',
        xclone_n_clones = None,
        set_spatial = False,
        xclone_plot = True,
        script_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/downsample_coverage/ds_50perc/cna_calling/xclone/scripts'
    )
