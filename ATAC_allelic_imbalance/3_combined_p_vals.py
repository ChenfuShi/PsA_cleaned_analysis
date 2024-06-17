#######
# this script takes the preprocessed combined individual p-value files and calculates grouped p-values for each variant

# note that we cannot provide raw genotypes for the patients, so this script is here only as a reference and it cannot be run
#######


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
import io
import functools
import sys
import pickle
from pandarallel import pandarallel
pandarallel.initialize(nb_workers = 8)

base_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_cleaned_analysis"

def calculate_p_vals(x):
    return scipy.stats.combine_pvalues(x)[1]

def calc_evidence(REF,ALT):
    tot_REF = sum(REF)
    tot_ALT = sum(ALT)
    ratio = tot_ALT/tot_REF
    tot_evidence = len(REF)
    return tot_REF,tot_ALT,ratio,tot_evidence

def apply_retrieve_p_vals(row,all_results):
    p_vals_greater = []
    p_vals_less = []
    ev_REF = []
    ev_ALT = []
    # Iterate through the dictionary of dataframes
    for key, df in all_results.items():
        # Check if the current dataframe contains a row with the same ID as the current sample
        match = df[(df['CHROM'] == row['CHROM']) & (df['POS'] == row['POS']) & 
            (df['ID'] == row['ID']) & (df['REF'] == row['REF']) & (df['ALT'] == row['ALT'])]
        if match.empty:
            continue
        # If there is a match, append the value of the specific column (key + "_A") to the list
        p_vals_greater.append(match[key + "_pval_greater"].values[0])
        p_vals_less.append(match[key + "_pval_less"].values[0])
        ev_REF.append(match[key + "_A"].values[0])
        ev_ALT.append(match[key + "_B"].values[0])
    # Apply the custom function to the list of values and add the result as a new column in the merged dataframe
    combined_p_val_greater = calculate_p_vals(p_vals_greater)
    combined_p_val_less = calculate_p_vals(p_vals_less)
    return combined_p_val_greater, combined_p_val_less, *calc_evidence(ev_REF,ev_ALT)

def identify_combined_p_vals_per_celltype(all_results):
    # Initialize an empty dataframe to store the merged data
    all_called_SNPs = pd.DataFrame()

    # Iterate through the dictionary of dataframes
    for key, df in all_results.items():
        # Select the first 5 columns from the current dataframe
        df = df.iloc[:, :5]
        # Append the current dataframe to the merged dataframe
        all_called_SNPs = all_called_SNPs.append(df, ignore_index=True)

    # Drop duplicates from the merged dataframe
    all_called_SNPs = all_called_SNPs.drop_duplicates()

    partially_applied_function = functools.partial(apply_retrieve_p_vals, all_results=all_results)

    all_called_SNPs[["combined_p_val_greater", "combined_p_val_less","tot_REF", "tot_ALT", "ratio", "n_pat"]] = all_called_SNPs.parallel_apply(partially_applied_function, axis = 1, result_type="expand")
    return all_called_SNPs

def main():
    metadata_ATAC = pd.read_csv(f"{base_dir}/metadata/cleaned_ATAC_metadata.csv", index_col=0)

    all_results = pickle.load(open(f"{base_dir}/ATAC_allelic_imbalance/all_p_val_called.pk", "rb"))
    all_results_CD4 = {k: all_results[k] for k in metadata_ATAC[metadata_ATAC["cell_type"] == "CD4"]["id"] if k in all_results}
    all_results_CD8 = {k: all_results[k] for k in metadata_ATAC[metadata_ATAC["cell_type"] == "CD8"]["id"] if k in all_results}
    all_results_CD4_all = {k: all_results[k] for k in metadata_ATAC[metadata_ATAC["cell_type"].isin(["CD4","CD4_SF"])]["id"] if k in all_results}
    all_results_CD8_all = {k: all_results[k] for k in metadata_ATAC[metadata_ATAC["cell_type"].isin(["CD8","CD8_SF"])]["id"] if k in all_results}
    all_results_CD4_SF = {k: all_results[k] for k in metadata_ATAC[metadata_ATAC["cell_type"].isin(["CD4_SF"])]["id"] if k in all_results}
    all_results_CD8_SF = {k: all_results[k] for k in metadata_ATAC[metadata_ATAC["cell_type"].isin(["CD8_SF"])]["id"] if k in all_results}

    index = int(sys.argv[1])

    # Define the list of functions to run
    functions_list = [identify_combined_p_vals_per_celltype(all_results_CD4),
                      identify_combined_p_vals_per_celltype(all_results_CD8),
                      identify_combined_p_vals_per_celltype(all_results_CD4_all),
                      identify_combined_p_vals_per_celltype(all_results_CD8_all),
                      identify_combined_p_vals_per_celltype(all_results),
                      identify_combined_p_vals_per_celltype(all_results_CD4_SF),
                      identify_combined_p_vals_per_celltype(all_results_CD8_SF),                      
                     ]
    # Run the function corresponding to the index
    result = functions_list[index]
    name = ["all_results_CD4","all_results_CD8","all_results_CD4_all","all_results_CD8_all","all_results","all_results_CD4_SF","all_results_CD8_SF"]
    # Save the result to a pickle file
    with open(f"{base_dir}/ATAC_allelic_imbalance/combined_p_vals_files/result_{name[index]}.pkl", "wb") as f:
        pickle.dump(result, f)




if __name__ == "__main__":
    main()





