######
# script that retrieves the counts for each loop from the two alleles for each sample from the bam file
# provided for reference only because we can't provide the bam files
######


import sys
import pandas as pd
import numpy as np
import os
import scipy
import pybedtools as pbed 
import subprocess as sub
from multiprocessing import Pool
import math
import statistics
import matplotlib.pyplot as plt
import seaborn as sns
os.makedirs("/mnt/iusers01/jw01/mdefscs4/scratch/temp_pybedtools/", exist_ok = True)
pbed.helpers.set_tempdir("/mnt/iusers01/jw01/mdefscs4/scratch/temp_pybedtools/")
bed_genome_file = "/mnt/iusers01/jw01/mdefscs4/hg38.genome"

plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['svg.fonttype'] = 'none'

def run_file(sample, patient, bam_A_file, bam_B_file):
    output_counts = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/calling_allele_specific/output_counts/" + sample + "_slop10kb.bedpe"
    bam_A = pbed.BedTool(bam_A_file).bam_to_bed(bedpe = True).sort()
    bam_B = pbed.BedTool(bam_B_file).bam_to_bed(bedpe = True).sort()

    interactions_f = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_combined_analysis/HiC_analysis/mustache_loops/merged_CD4_CD8_loops_reformat.bed"
    interactions = pd.read_csv(interactions_f, sep = "\t", header = None)
    res_A = bam_A.pair_to_pair(b = pbed.BedTool(interactions_f), type = "both", slop = 10000, **{"is": True})
    A_df = res_A.to_dataframe(header = None, disable_auto_names=True)
    A_df = A_df.groupby([10,11,12,13,14,15,16,17]).count().reset_index()[[10,11,12,13,14,15,16,17,0]]
    A_df = A_df.rename({0: "counts_A"}, axis = 1)
    res_B = bam_B.pair_to_pair(b = pbed.BedTool(interactions_f), type = "both", slop = 10000, **{"is": True})
    B_df = res_B.to_dataframe(header = None, disable_auto_names=True)
    B_df = B_df.groupby([10,11,12,13,14,15,16,17]).count().reset_index()[[10,11,12,13,14,15,16,17,0]]
    B_df = B_df.rename({0: "counts_B"}, axis = 1)
    interactions.merge(A_df,
    left_on = [0,1,2,3,4,5,6,7], 
    right_on = [10,11,12,13,14,15,16,17], 
    how = "outer")[[0,1,2,3,4,5,6,7,"counts_A"]].merge(B_df,
    left_on = [0,1,2,3,4,5,6,7], 
    right_on = [10,11,12,13,14,15,16,17], 
    how = "outer")[[0,1,2,3,4,5,6,7,"counts_A","counts_B"]].fillna(0).to_csv(output_counts,sep="\t",header = False, index = False)
    



# Get command line arguments
sample = sys.argv[1]
patient = sys.argv[2]

bam_A_file = f"/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/allele_splitting/splitted_bams/{sample}_G1_reads_byname.bam"
bam_B_file = f"/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/allele_splitting/splitted_bams/{sample}_G2_reads_byname.bam"


run_file(sample, patient, bam_A_file, bam_B_file)