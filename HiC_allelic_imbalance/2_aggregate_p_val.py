######
# script that takes the counts from the previous step and generates p-values for each loop
# provided for reference only because we can't provide the bam files
######


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
import statsmodels.stats.multitest
import gzip
import subprocess
import vcf
import functools
os.makedirs("/mnt/iusers01/jw01/mdefscs4/scratch/temp_pybedtools/", exist_ok = True)
pbed.helpers.set_tempdir("/mnt/iusers01/jw01/mdefscs4/scratch/temp_pybedtools/")
bed_genome_file = "/mnt/iusers01/jw01/mdefscs4/hg38.genome"

plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['svg.fonttype'] = 'none'

base_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_cleaned_analysis"

metadata_hic = pd.read_csv(f"{base_dir}/metadata/cleaned_HiC_metadata.csv", index_col=0)

def retrieve_datatable(samples):
    sample = samples[0]
    counts = pd.read_csv(f"/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/calling_allele_specific/output_counts/{sample}_slop10kb.bedpe", sep ="\t", header = None)
    counts = counts[[0,1,2,3,4,5,6,7]]
    counts[6] = counts.index
    for sample in samples:
        counts_B = pd.read_csv(f"/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/calling_allele_specific/output_counts/{sample}_slop10kb.bedpe", sep ="\t", header = None)
        counts_B = counts_B.rename(columns={8:sample + "_A",9:sample + "_B"})
        counts_B[6] = counts_B.index
        counts = counts.merge(counts_B)
    counts_pbed = pbed.BedTool.from_dataframe(counts)
    return counts, counts_pbed


samples_CD4 = metadata_hic[metadata_hic["cell_type"] == "CD4"]["folder_name"].to_list()
samples_CD8 = metadata_hic[metadata_hic["cell_type"] == "CD8"]["folder_name"].to_list()
counts_CD4, counts_pbed_CD4 = retrieve_datatable(samples_CD4)
counts_CD8, counts_pbed_CD8 = retrieve_datatable(samples_CD8)

def get_if_good(r):
    tot = 0
    for s in r.samples:
        if s.data.GT in ["0|1","1|0"]:
            tot = tot + 1
    if tot >= 2:
        return True
    else:
        return False

def extract_id(string):
    return os.path.split(string)[1].split("_")[0]

def run_calling(chrom, counts, counts_pbed, samples):

    os.makedirs(f"/mnt/iusers01/jw01/mdefscs4/scratch/temp_pybedtools_{chrom}/", exist_ok = True)
    pbed.helpers.set_tempdir(f"/mnt/iusers01/jw01/mdefscs4/scratch/temp_pybedtools_{chrom}/")

    valid_patients = metadata_hic[metadata_hic["folder_name"].isin(samples)]["patient"].to_list()
    vcf_f = pbed.BedTool(f"/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/genotype_calling/called_genotypes/HAPCUT_phased_merged_chr{chrom}_annotated_filtered.vcf.gz")
    loop_anchors = pbed.BedTool("/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/calling_allele_specific/loop_anchors_merged.bed")
    filtered_vcf = vcf_f.intersect(loop_anchors, u = True, header = True)

    variants = vcf.Reader(filename=filtered_vcf.fn)
    counts_pbed_chrom = counts_pbed.intersect(pbed.BedTool(f"chr{chrom}\t1\t300000000", from_string = True), u = True)

    dataset = {}
    for record in variants:
        if get_if_good(record):
            variant_pos = pbed.BedTool(f"{record.CHROM}\t{record.POS - 1}\t{record.POS}", from_string = True)
            useful_interactions = counts_pbed_chrom.pairtobed(variant_pos).to_dataframe(header = None, disable_auto_names=True).iloc[:,:-3]
            if len(useful_interactions) < 1:
                continue
            useful_interactions.columns = counts.columns
            for id, line in useful_interactions.iterrows():
                ref_idx = f"{record.ID}_{line[6]}"
                idx_data = {}
                for s in record.samples:
                    ID_pat = extract_id(str(s.sample)) 
                    if ID_pat in valid_patients:
                        ID_sample = metadata_hic[(metadata_hic["patient"] == ID_pat) & (metadata_hic["folder_name"].isin(samples))]["folder_name"].to_list()[0]
                        if line[f"{ID_sample}_A"]+line[f"{ID_sample}_B"] > 10 and line[f"{ID_sample}_A"] > 1 and line[f"{ID_sample}_B"] > 1:
                            if s.data.GT == "1|0":
                                idx_data[f"{ID_sample}_1"] = line[f"{ID_sample}_A"]
                                idx_data[f"{ID_sample}_0"] = line[f"{ID_sample}_B"]
                                idx_data[f"{ID_sample}_pval_greater"] = scipy.stats.binom_test(idx_data[f"{ID_sample}_1"],idx_data[f"{ID_sample}_1"] + idx_data[f"{ID_sample}_0"],0.5,"greater")
                                idx_data[f"{ID_sample}_pval_less"] = scipy.stats.binom_test(idx_data[f"{ID_sample}_1"],idx_data[f"{ID_sample}_1"] + idx_data[f"{ID_sample}_0"],0.5,"less")
                            elif s.data.GT == "0|1":
                                idx_data[f"{ID_sample}_1"] = line[f"{ID_sample}_B"]
                                idx_data[f"{ID_sample}_0"] = line[f"{ID_sample}_A"]
                                idx_data[f"{ID_sample}_pval_greater"] = scipy.stats.binom_test(idx_data[f"{ID_sample}_1"],idx_data[f"{ID_sample}_1"] + idx_data[f"{ID_sample}_0"],0.5,"greater")
                                idx_data[f"{ID_sample}_pval_less"] = scipy.stats.binom_test(idx_data[f"{ID_sample}_1"],idx_data[f"{ID_sample}_1"] + idx_data[f"{ID_sample}_0"],0.5,"less")
                if idx_data:
                    keys_pval = [x for x in idx_data.keys() if "pval_greater" in x]
                    idx_data["combined_p_val_greater"] = scipy.stats.combine_pvalues([idx_data[x] for x in keys_pval])[1]
                    keys_pval = [x for x in idx_data.keys() if "pval_less" in x]
                    idx_data["combined_p_val_less"] = scipy.stats.combine_pvalues([idx_data[x] for x in keys_pval])[1]
                    dataset[ref_idx] = idx_data
    subprocess.run(f"find /mnt/iusers01/jw01/mdefscs4/scratch/temp_pybedtools_{chrom}/ -mmin +59 -type f -exec rm -f {'{}'} \;", shell = True)
    return dataset


partial_run_calling_CD4 = functools.partial(run_calling, counts=counts_CD4, counts_pbed=counts_pbed_CD4, samples=samples_CD4)
partial_run_calling_CD8 = functools.partial(run_calling, counts=counts_CD8, counts_pbed=counts_pbed_CD8, samples=samples_CD8)

import sys
sample = sys.argv[1]
import pickle

if sample == "CD4":
    with Pool(18) as p:
        data_CD4 = p.map(partial_run_calling_CD4, list(range(1,23)))
    pickle.dump(data_CD4, open("/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/calling_allele_specific/aggregated_data_CD4_slop10kb_separatepval.pk", "wb"))

elif sample == "CD8":
    with Pool(18) as p:
        data_CD8 = p.map(partial_run_calling_CD8, list(range(1,23)))
    pickle.dump(data_CD8, open("/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/calling_allele_specific/aggregated_data_CD8_slop10kb_separatepval.pk", "wb"))
