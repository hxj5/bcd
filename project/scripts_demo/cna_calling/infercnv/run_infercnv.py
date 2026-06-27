# run_infercnv.py


import os
import stat
import subprocess
import sys
import time
from logging import error


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



if __name__ == "__main__":
    run_infercnv(
        sample_id = 'HCC3N_simu',
        matrix_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/downsample_coverage/ds_50perc/gen_data/rdr_starsolo',
        out_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/downsample_coverage/ds_50perc/cna_calling/infercnv',
        cell_anno_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/base/data/matrix/spot_anno.tsv',
        ref_cell_types = 'normal',
        gene_anno_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/refapp/infercnv/infercnv.genes.tsv',
        ncores = 10,
        cutoff = 0.1,
        png_dpi = 300,
        script_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/downsample_coverage/ds_50perc/cna_calling/infercnv/scripts/run_infercnv.R'
    )
