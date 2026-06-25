# gen_features.py

import pandas as pd
import sys, os
import scReadSim.Utility as Utility
import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount
import scReadSim.scRNA_GenerateBAM as scRNA_GenerateBAM


root_dir = "/groups/cgsd/xianjie/projects/cna-benchmark/HCC3N_600spot/simulator_evaluation/base/gen_data/scdesign2_screadsim"
out_dir = os.path.join(root_dir, "simu/matrix")


INPUT_cells_barcode_file = os.path.join(root_dir, 'data/barcodes.tsv')
INPUT_bamfile = os.path.join(root_dir, 'simu/pp/processed.bam')
INPUT_genome_size_file = os.path.join(root_dir, 'data/hg38_chrom_sizes.tsv')
INPUT_genome_annotation = os.path.join(root_dir, "data/genes.gtf")

samtools_directory = "/home/xianjie/.anaconda3/envs/SCSC/bin"
macs3_directory = "/home/xianjie/.anaconda3/envs/SCSC/bin/macs3"
bedtools_directory = "/home/xianjie/.anaconda3/envs/SCSC/bin"
seqtk_directory = "/home/xianjie/.anaconda3/envs/SCSC/bin/seqtk"
fgbio_jarfile = "/home/xianjie/.anaconda3/envs/SCSC/bin/fgbio"


# Generate features
outdirectory = out_dir
os.makedirs(outdirectory, exist_ok = True)

Utility.scRNA_CreateFeatureSets(
    INPUT_bamfile = INPUT_bamfile,
    samtools_directory = samtools_directory,
    bedtools_directory = bedtools_directory,
    outdirectory = outdirectory,
    genome_annotation = INPUT_genome_annotation,
    genome_size_file = INPUT_genome_size_file
)

# Keep only chr1-22.
fn = os.path.join(out_dir, 'scReadSim.Gene.bed')
out_fn = os.path.join(out_dir, 'scReadSim.Gene.chr1-22.bed')
target_chroms = [str(i) for i in range(1, 23)]
df = pd.read_csv(fn, header = None, sep = '\t')
df.columns = ['chrom', 'start', 'end']
df['chrom'] = df['chrom'].astype(str)
df = df.loc[df['chrom'].isin(target_chroms)].copy()
df.to_csv(out_fn, header = False, sep = '\t', index = False)

fn = os.path.join(out_dir, 'scReadSim.InterGene.bed')
out_fn = os.path.join(out_dir, 'scReadSim.InterGene.chr1-22.bed')
target_chroms = [str(i) for i in range(1, 23)]
df = pd.read_csv(fn, header = None, sep = '\t')
df.columns = ['chrom', 'start', 'end']
df['chrom'] = df['chrom'].astype(str)
df = df.loc[df['chrom'].isin(target_chroms)].copy()
df.to_csv(out_fn, header = False, sep = '\t', index = False)

print("[CreateFeatureSets] All Done!")
