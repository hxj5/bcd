# run_calicost.py



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



def run_calicost(
    sample_id,
    out_dir,
    ncores,
    genome,
    eagle_fn,
    snp_vcf_fn,
    panel_dir,
    calicost_dir,
    bam_fn = None,
    spaceranger_dir = None,
    bamlist = None,
    normalidx_file = None,
    script_dir = None
):
    """
    Run CalicoST for CNA calling, skipping `purity` step:
    1. Preprocessing (genotyping and phasing)
    2. Run CalicoST main (CNA calling)

    Supports both single sample and multiple samples.

    Parameters
    ----------
    sample_id : str or list
        Sample ID(s). For single sample, provide a string. For multiple samples,
        use `bamlist`.
    out_dir : str
        Output directory.
    ncores : int
        Number of cores.
    genome : str
        Genome version (e.g., "hg38", "hg19").
    eagle_fn : str
        Path to Eagle executable.
    snp_vcf_fn : str
        Path to SNP VCF file.
    panel_dir : str
        Path to phasing panel directory.
    calicost_dir : str
        Path to CalicoST installation directory.
    bam_fn : str
        Path to BAM file. For single sample, provide a string. For multiple
        samples, use `bamlist`.
    spaceranger_dir : str
        Path to spaceranger output directory (for single sample).
    bamlist : list of tuples, optional
        For multiple samples, list of (bam_path, sample_id, spaceranger_dir).
        If provided, overrides bam_fn and sample_id.
    normalidx_file : str, optional
        The path to the file containing the indices of normal spots in the 
        spatial transcriptomics data. 
        Each line is a single index without header.
    script_dir : str, optional
        Directory for scripts. If None, will be set to out_dir/scripts.

    Returns
    -------
    int
        Return code from the command execution.
    """
    # default settings.
    skip_preprocessing = False


    # init.
    assert genome in ["hg19", "hg38"]
    genome_GRCh = "GRCh38" if genome == 'hg38' else "GRCh37"
    geneticmap_file = os.path.join(
        calicost_dir, 
        "%s_resources" % genome_GRCh, 
        "genetic_map_%s_merged.tab.gz" % genome_GRCh
    )
    hgtable_file = os.path.join(
        calicost_dir, 
        "%s_resources" % genome_GRCh,
        "hgTables_%s_gencode.txt" % genome
    )
    filtergenelist_file = os.path.join(
        calicost_dir, 
        "%s_resources" % genome_GRCh,
        "ig_gene_list.txt"
    )
    filterregion_file = os.path.join(
        calicost_dir, 
        "%s_resources" % genome_GRCh,
        "HLA_regions.bed"
    )
    
    os.makedirs(out_dir, exist_ok = True)

    if script_dir is None:
        script_dir = os.path.join(out_dir, "scripts")
    os.makedirs(script_dir, exist_ok = True)

    config_dir = os.path.join(out_dir, "config")
    os.makedirs(config_dir, exist_ok = True)

    snp_dir = os.path.join(out_dir, "snp")
    os.makedirs(snp_dir, exist_ok = True)

    is_multi, bamlist_entries, single_spaceranger_dir = None, None, None
    if bamlist is None:
        is_multi = False
        assert os.path.exists(spaceranger_dir)
        bamlist_entries = [(bam_fn, sample_id, spaceranger_dir)]
        single_spaceranger_dir = spaceranger_dir
    else:
        is_multi = True
        bamlist_entries = bamlist
        single_spaceranger_dir = bamlist_entries[0][2]

    bamlist_fn = os.path.join(config_dir, "bamlist.tsv")
    with open(bamlist_fn, "w") as fp:
        for bam_path, sid, sp_dir in bamlist_entries:
            fp.write("%s\t%s\t%s\n" % (bam_path, sid, sp_dir))


    # Step 1: Preprocessing (genotyping and phasing)
    print("[D] pp ...")
    if not skip_preprocessing:
        config_pp_fn = os.path.join(config_dir, "config.yaml")
        with open(config_pp_fn, "w") as fp:
            fp.write("# path to executables or their parent directories\n")
            fp.write("calicost_dir: %s\n" % calicost_dir)
            fp.write("eagledir: %s\n" % os.path.dirname(eagle_fn))
            fp.write("\n")
            fp.write("# running parameters\n")
            fp.write("# samtools sort (only used when joingly calling from multiple slices)\n")
            fp.write("samtools_sorting_mem: \"4G\"\n")
            fp.write("# cellsnp-lite\n")
            fp.write("UMItag: \"Auto\"\n")
            fp.write("cellTAG: \"CB\"\n")
            fp.write("nthreads_cellsnplite: %d\n" % ncores)
            fp.write("region_vcf: %s\n" % snp_vcf_fn)
            fp.write("# Eagle phasing\n")
            fp.write("phasing_panel: %s\n" % panel_dir)
            fp.write("chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]\n")
            fp.write("\n")
            fp.write("# input\n")
            if is_multi:
                fp.write("bamlist: %s\n" % bamlist_fn)
            else:
                fp.write("spaceranger_dir: %s\n" % single_spaceranger_dir)
            fp.write("\n")
            fp.write("# output\n")
            fp.write("output_snpinfo: %s\n" % snp_dir)

        # Create script for preprocessing
        script_pp_fn = os.path.join(script_dir, "calicost.pp.sh")
        s  = '''#!/bin/bash\n'''
        s += '''# run CalicoST preprocessing.\n'''
        s += '''\n'''
        s += '''set -eux'''
        s += '''\n'''
        s += '''snakemake        \\\n'''
        s += '''    --cores  %d  \\\n''' % ncores
        s += '''    --configfile  %s  \\\n''' % config_pp_fn
        s += '''    --snakefile   %s  all\n''' % os.path.join(calicost_dir, "calicost.smk")
        s += '''\n'''
        s += '''echo "[CalicoST Preprocessing] All Done!"\n'''
        s += '''\n'''

        with open(script_pp_fn, "w") as fp:
            fp.write(s)
        set_file_exe(script_pp_fn)

        cmd = "/usr/bin/time -v %s" % script_pp_fn
        ret, outs, errs = run_cmd(cmd)
        if ret != 0:
            return(ret)


    # Step 2: Run CalicoST main (CNA calling)
    print("[D] cna ...")
    config_cna_fn = os.path.join(config_dir, "configuration_cna")
    with open(config_cna_fn, "w") as fp:
        fp.write("\n")
        if is_multi:
            fp.write("input_filelist : %s\n" % bamlist_fn)
        else:
            fp.write("spaceranger_dir : %s\n" % single_spaceranger_dir)
        fp.write("snp_dir : %s\n" % snp_dir)
        fp.write("output_dir : %s\n" % out_dir)
        fp.write("\n")
        fp.write("# supporting files and preprocessing arguments\n")
        fp.write("geneticmap_file : %s\n" % geneticmap_file)
        fp.write("hgtable_file : %s\n" % hgtable_file)
        if normalidx_file is None:
            fp.write("normalidx_file : None\n")
        else:
            fp.write("normalidx_file : %s\n" % normalidx_file)
        fp.write("tumorprop_file : None\n")
        fp.write("alignment_files : \n")
        fp.write("supervision_clone_file : None\n")
        fp.write("filtergenelist_file : %s\n" % filtergenelist_file)
        fp.write("filterregion_file : %s\n" % filterregion_file)
        fp.write("secondary_min_umi : 300\n")
        fp.write("bafonly : False\n")
        fp.write("\n")
        fp.write("# phase switch probability\n")
        fp.write("nu : 1.0\n")
        fp.write("logphase_shift : -2.0\n")
        fp.write("npart_phasing : 3\n")
        fp.write("\n")
        fp.write("# HMRF configurations\n")
        fp.write("n_clones : 3\n")
        fp.write("n_clones_rdr : 2\n")
        fp.write("min_spots_per_clone : 100\n")
        fp.write("min_avgumi_per_clone : 10\n")
        fp.write("maxspots_pooling : 7\n")
        fp.write("tumorprop_threshold : 0.5\n")
        fp.write("max_iter_outer : 20\n")
        fp.write("nodepotential : weighted_sum\n")
        fp.write("initialization_method : rectangle\n")
        fp.write("num_hmrf_initialization_start : 0\n")
        fp.write("num_hmrf_initialization_end : 1\n")
        fp.write("spatial_weight : 1.0\n")
        fp.write("construct_adjacency_method : hexagon\n")
        fp.write("construct_adjacency_w : 1.0\n")
        fp.write("\n")
        fp.write("# HMM configurations\n")
        fp.write("n_states : 7\n")
        fp.write("params : smp\n")
        fp.write("t : 1-1e-4\n")
        fp.write("t_phaseing : 0.9999\n")
        fp.write("fix_NB_dispersion : False\n")
        fp.write("shared_NB_dispersion : True\n")
        fp.write("fix_BB_dispersion : False\n")
        fp.write("shared_BB_dispersion : True\n")
        fp.write("max_iter : 30\n")
        fp.write("tol : 0.0001\n")
        fp.write("gmm_random_state : 0\n")
        fp.write("np_threshold : 1.0\n")
        fp.write("np_eventminlen : 10\n")
        fp.write("\n")
        fp.write("# integer copy number\n")
        fp.write("nonbalance_bafdist : 1.0\n")
        fp.write("nondiploid_rdrdist : 10.0\n")

    # Create Python script to run CalicoST
    script_fn = os.path.join(script_dir, "calicost.call.sh")
    s  = '''#!/bin/bash\n'''
    s += '''# run CalicoST CNA calling.\n'''
    s += '''\n'''
    s += '''set -eux'''
    s += '''\n'''
    s += '''OMP_NUM_THREADS=1\n'''
    s += '''\n'''
    s += '''python %s  \\\n''' % os.path.join(calicost_dir, "src/calicost/calicost_main.py")
    s += '''    -c  %s\n''' % config_cna_fn
    s += '''\n'''
    s += '''echo "[CalicoST CNA] All Done!"\n'''
    s += '''\n'''

    with open(script_fn, "w") as fp:
        fp.write(s)
    set_file_exe(script_fn)

    cmd = "/usr/bin/time -v %s" % script_fn
    ret, outs, errs = run_cmd(cmd)
    print("[D] end ...")
    return(ret)



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
    run_calicost(
        sample_id = 'HCC3N_simu',
        out_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1r1n1t/vary_tumor_prop/700n300t/cna_calling/calicost_skippurity',
        ncores = 10,
        genome = 'hg38',
        eagle_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/refapp/eagle/Eagle_v2.4.1/eagle',
        snp_vcf_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/refapp/phasing/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz',
        panel_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/refapp/phasing/1000G_hg38',
        calicost_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/base/data/refapp/calicost/CalicoST',
        bam_fn = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1r1n1t/vary_tumor_prop/700n300t/gen_data/ps_spaceranger/possorted_genome_bam.bam',
        spaceranger_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1r1n1t/vary_tumor_prop/700n300t/gen_data/ps_spaceranger',
        bamlist = None,
        normalidx_file = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1r1n1t/vary_tumor_prop/700n300t/gen_data/ps_spaceranger/normalidx_file.tsv',
        script_dir = '/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/cna_prediction/1r1n1t/vary_tumor_prop/700n300t/cna_calling/calicost_skippurity/scripts'
    )
