import pandas as pd
import numpy as np
import hicstraw 
from multiprocessing import Pool
from functools import partial
import glob
import os
import plotly.express as px
import math
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from pandarallel import pandarallel
import cooler
import cooltools
import pybedtools as pbed
pandarallel.initialize()
from scipy import stats, special
from statsmodels.stats import multitest
import statsmodels.api as sm
import statsmodels.formula.api as smf
import plotly.io as pio
import seaborn as sns
import numba as nb
import bioframe
os.makedirs("/mnt/iusers01/jw01/mdefscs4/scratch/temp_pybedtools/", exist_ok = True)
pbed.helpers.set_tempdir("/mnt/iusers01/jw01/mdefscs4/scratch/temp_pybedtools/")
bed_genome_file = "/mnt/iusers01/jw01/mdefscs4/hg38.genome"

plt.rcParams['svg.fonttype'] = 'none'



def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols
        
def retrieve_df_correction_factors(sample, chrom, cvd_base): 
    """ override original to load the precomputed one"""
    sample_correction = pd.read_pickle(f"/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_combined_analysis/HiC_analysis/temp_output/correction_factors/{sample}.pk")
    if chrom:
        sample_correction = sample_correction[sample_correction["region1"] == chrom]
    return sample_correction

def normalize_sample(mat, sample_correction, resolution = 5000):
    for i in range(mat.shape[0]):
        idx = kth_diag_indices(mat, k = i)
        idx_2 = kth_diag_indices(mat, k = -i)
        if resolution == 2500:
            mat[idx] = mat[idx] * sample_correction[sample_correction["dist"] == math.floor(i/4) + 1]["correction_factor"].values
            mat[idx_2] = mat[idx_2] * sample_correction[sample_correction["dist"] == math.floor(i/4) + 1]["correction_factor"].values
        elif resolution == 5000:
            mat[idx] = mat[idx] * sample_correction[sample_correction["dist"] == math.floor(i/2) + 1]["correction_factor"].values
            mat[idx_2] = mat[idx_2] * sample_correction[sample_correction["dist"] == math.floor(i/2) + 1]["correction_factor"].values
        elif resolution == 10000:
            mat[idx] = mat[idx] * sample_correction[sample_correction["dist"] == i]["correction_factor"].values
            mat[idx_2] = mat[idx_2] * sample_correction[sample_correction["dist"] == i]["correction_factor"].values    
        elif resolution == 25000:       
            mat[idx] = mat[idx] * sample_correction[sample_correction["dist"] == math.floor(i*2.5)]["correction_factor"].values
            mat[idx_2] = mat[idx_2] * sample_correction[sample_correction["dist"] == math.floor(i*2.5)]["correction_factor"].values
        elif resolution == 50000:       
            mat[idx] = mat[idx] * sample_correction[sample_correction["dist"] == math.floor(i*5)]["correction_factor"].values
            mat[idx_2] = mat[idx_2] * sample_correction[sample_correction["dist"] == math.floor(i*5)]["correction_factor"].values
    return mat

def get_data_sample(sample,chrom,region_start,region_end,resolution = 5000):
    hic_file = f"/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_HiC_analyzer/juiceboxes/{sample}/{sample}.allValidPairs.hic"
    hic = hicstraw.HiCFile(hic_file)
    mzd = hic.getMatrixZoomData(chrom, chrom, "observed", "SCALE", "BP", resolution)
    data = mzd.getRecordsAsMatrix(region_start, region_end, region_start, region_end)
    sample_correction = retrieve_df_correction_factors(sample, chrom, cvd_base = None) 
    data = normalize_sample(data, sample_correction)    
    return data

def get_region_for_all(required_samples,chrom,region_start,region_end,resolution = 5000):
    pool = Pool(28)
    get_data_sample_part = partial(get_data_sample, chrom=chrom,region_start=region_start,region_end=region_end,resolution=resolution)
    all_data = list(pool.map(get_data_sample_part, required_samples))
    orig_shape = all_data[0].shape
    matrix_full = np.array(all_data).sum(axis = 0)
    all_data = [x.flatten() for x in all_data]
    return all_data, orig_shape, matrix_full

def retrieve_stats_scipy(X, y):
    y = np.array(y)
    results = np.empty((5, X.shape[1]), dtype = np.float64)
    for i in range(X.shape[1]):
        if np.all(X[:,i] == X[:,i][0]):
            results[:,i] = 0,0,0,1,0
        else:
            results[:,i] = stats.linregress(X[:,i], y) # slope, intercept, r_value, p_value, std_err
    return results

def retrieve_stats_statsmodels(X, y):
    y = np.array(y)
    results = np.empty((5, X.shape[1]), dtype = np.float64)
    for i in range(X.shape[1]):
        if np.all(X[:,i] == X[:,i][0]):
            results[:,i] = 0,0,0,0,1
        else:
            loc_X = X[:,i]
            loc_X = sm.add_constant(loc_X)
            model = sm.OLS(y, loc_X).fit()
            results[:,i] = model.params[1], model.params[0], model.rsquared, model.rsquared_adj, model.f_pvalue
    return results

def myfloor(x, base=2500):
    return base * math.floor(x/base)
def myceil(x, base=2500):
    return base * math.ceil(x/base)


def add_numbers(df,local_genes):
    def _assign_number(group):
        # get the gene name
        gene = group.iloc[0]['gene_name']
        # create a dictionary that maps each gene to a unique number
        shape = gene_list.shape
        seq = np.tile(np.arange(0.5, 3.5, 1.4), int(np.prod(shape)/3) + 1)[:np.prod(shape)]*0.82
        gene_map = {gene_list[i]: k for i,k in enumerate(seq)}
        return pd.Series({'position': gene_map[gene]})
    gene_list = np.unique(df['gene_name'].to_list() + local_genes['gene_name'].to_list())
    df_grouped = df.groupby('gene_name').apply(_assign_number)
    df_merged = pd.merge(df, df_grouped, on='gene_name')
    df_grouped = local_genes.groupby('gene_name').apply(_assign_number)
    local_genes = pd.merge(local_genes, df_grouped, on='gene_name')
    return df_merged, local_genes
