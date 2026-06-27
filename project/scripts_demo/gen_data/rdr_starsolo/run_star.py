# run_star.py


import numpy as np
import os
import pandas as pd
import shutil
import stat
import subprocess
import sys



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



if __name__ == "__main__":
    run_starsolo(
        bam_fn = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/downsample_coverage/ds_50perc/gen_data/subset_bam/possorted.bam",
        genome_dir = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/refseq/star/genome_index",
        cb_whitelist_fn = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/base/data/starsolo.barcodes.pseudo_whitelist.tsv",
        out_dir = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/downsample_coverage/ds_50perc/gen_data/rdr_starsolo",
        ncores = 10,
        read_file_type = 'SAM SE',
        read_file_command = 'samtools view -F 0x100',
        solo_strand = 'Forward',
        solo_type = 'CB_UMI_Simple',
        solo_barcode_tags = 'CR UR',
        solo_qual_tags = 'CY UY',
        solo_umi_len = 12,
        out_sam_type = 'BAM Unsorted',
        script_fn = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1n1t_s1600/downsample_coverage/ds_50perc/gen_data/rdr_starsolo/scripts/run_star.sh"
    )

