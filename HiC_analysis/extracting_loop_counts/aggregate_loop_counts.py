import pandas as pd
import numpy as np

import glob
import os

import cooler
import cooltools

##########################
# script for the aggregation in batch of all the precomputed loop counts
# takes the counts and the cool files and calculates how to normalize the loop counts by the distance distribution


##########################

base_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_cleaned_analysis"
metadata_HiC_file = f"{base_dir}/metadata/cleaned_HiC_metadata.csv"

bias_res = 10000
reference_file = "PSA4920_CD4_ARIMA" # this can be any file
loop_count_dir = f"{base_dir}/HiC_analysis/extracting_loop_counts/counts/"

def get_distance(row):
    """calculates the distance bin of this particular loop so that it can be corrected by the decay factor"""
    return round(abs(row["A_start"] - row["B_start"])/bias_res)

def retrieve_df_correction_factors(sample, chrom = None, cvd_base = None): 
    """ load precomputed correction factors"""
    sample_correction = pd.read_pickle(f"{base_dir}/HiC_analysis/extracting_loop_counts/correction_factors/{sample}.pk")
    if chrom:
        sample_correction = sample_correction[sample_correction["region1"] == chrom]
    return sample_correction


if __name__=="__main__":

    metadata_hic = pd.read_csv(metadata_HiC_file, index_col = 0)

    files_CD8_8_SF = glob.glob(loop_count_dir + "*CD8_CD8SF")
    loops_analysed = pd.read_csv(files_CD8_8_SF[0], sep = "\t", dtype = {"chrA":str,"chrB":str})[['chrA', 'A_start', 'A_end', 'chrB', 'B_start', 'B_end', 'FDR', 'DETECTION_SCALE']]
    loops_analysed["distance_bin"] = loops_analysed.apply(get_distance, axis = 1)
    for sample in files_CD8_8_SF:
        basename = os.path.basename(sample).split("_counts_")[0]
        sample_correction = retrieve_df_correction_factors(basename)
        df = pd.read_csv(sample, sep = "\t", dtype = {"chrA":str,"chrB":str})
        df["distance_bin"] = df.apply(get_distance, axis = 1)
        df = df.merge(sample_correction, left_on = ["chrA","distance_bin"], right_on = ["region1", "dist"], how = "left")
        df["extracted_counts"] = df["extracted_counts"] * df["correction_factor"]
        loops_analysed[basename] = df[["extracted_counts"]]

    loops_analysed.to_pickle(f"{base_dir}/HiC_analysis/extracting_loop_counts/aggregated_counts/aggregated_normalized_loops_CD8_CD8SF.pk")


    files_CD4_8 = glob.glob(loop_count_dir + "*CD8_CD4")
    loops_analysed = pd.read_csv(files_CD4_8[0], sep = "\t", dtype = {"chrA":str,"chrB":str})[['chrA', 'A_start', 'A_end', 'chrB', 'B_start', 'B_end', 'FDR', 'DETECTION_SCALE']]
    loops_analysed["distance_bin"] = loops_analysed.apply(get_distance, axis = 1)
    for sample in files_CD4_8:
        basename = os.path.basename(sample).split("_counts_")[0]
        sample_correction = retrieve_df_correction_factors(basename)
        df = pd.read_csv(sample, sep = "\t", dtype = {"chrA":str,"chrB":str})
        df["distance_bin"] = df.apply(get_distance, axis = 1)
        df = df.merge(sample_correction, left_on = ["chrA","distance_bin"], right_on = ["region1", "dist"], how = "left")
        df["extracted_counts"] = df["extracted_counts"] * df["correction_factor"]
        loops_analysed[basename] = df[["extracted_counts"]]

    loops_analysed.to_pickle(f"{base_dir}/HiC_analysis/extracting_loop_counts/aggregated_counts/aggregated_normalized_loops_CD4_CD8.pk")

