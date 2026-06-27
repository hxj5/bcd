# run_copykat.py


import os
import stat
import subprocess
import sys
import time
from logging import error


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
    run_copykat(
        sample_id = 'HCC3N_simu',
        matrix_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1r1n1t/vary_tumor_prop/700n300t/gen_data/matrix_noref',
        out_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1r1n1t/vary_tumor_prop/700n300t/cna_calling/copykat_noref',
        cell_anno_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1r1n1t/vary_tumor_prop/700n300t/gen_data/matrix_noref/spot_anno.tsv',
        ref_cell_types = None,
        ncores = 10,
        script_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1r1n1t/vary_tumor_prop/700n300t/cna_calling/copykat_noref/scripts/run_copykat.R'
    )
