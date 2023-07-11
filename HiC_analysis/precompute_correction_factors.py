##############################
# This script needs to be run to pre-compute the cross-sample correction factors for the hic
# This is needed to really increase the speed of the computation of the plots

##############################
import pandas as pd
import numpy as np
import cooltools
import cooler
import bioframe

base_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_cleaned_analysis"
metadata_HiC_file = f"{base_dir}/metadata/cleaned_HiC_metadata.csv"

def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols
        
def retrieve_expected_cis(sample, chrom = None):
    clr = cooler.Cooler(f'/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_HiC_analyzer/cool/{sample}.mcool::/resolutions/'+str(10000))
    if chrom:
        t = pd.DataFrame({"chrom":[chrom], "start":[0], "end":[bioframe.fetch_chromsizes('hg38')[f"chr{chrom}"]], "name":[chrom]})
    else:
        t = None
    cvd = cooltools.expected_cis(
        clr=clr,
        aggregate_smoothed=False,
        nproc=16, #if you do not have multiple cores available, set to 1
        clr_weight_name=None,  
        ignore_diags = 1,
        view_df=t,
        intra_only=True,
        smooth = True,
        )
    cvd = cvd[cvd["dist"] <= 200]
    return cvd

bias_res = 10000
reference_file = "PSA4920_CD4_ARIMA"
def run_precompute():
    def _retrieve_df_correction_factors(sample, chrom, cvd_base):
        cvd_sample = retrieve_expected_cis(sample, chrom)
        cvd_sample = cvd_sample.merge(cvd_base[["region1","dist","count.avg.smoothed"]] , on = ["region1","dist"], validate = "one_to_one")
        cvd_sample["correction_factor"] = cvd_sample["count.avg.smoothed_y"] / cvd_sample["count.avg.smoothed_x"]
        return cvd_sample[["region1", "dist", "correction_factor"]]
    metadata_hic = pd.read_csv(metadata_HiC_file, index_col = 0)

    cvd_base = retrieve_expected_cis(reference_file)
    # precompute sample corrections
    for sample in metadata_hic["folder_name"].to_list():
        sample_correction = _retrieve_df_correction_factors(sample, None, cvd_base)
        sample_correction.to_pickle(f"{base_dir}/HiC_analysis/extracting_loop_counts/correction_factors/{sample}.pk")

    # we also want to run the precompute for the cell lines and the merged files
    extra_files = ["all_jurkat","MyLa_all_old_ARIMA","big_merger_all_CD8_SF_oct22","big_merger_all_CD8_oct22","big_merger_all_CD4_oct22"]
    for sample in extra_files:
        sample_correction = _retrieve_df_correction_factors(sample, None, cvd_base)
        sample_correction.to_pickle(f"{base_dir}/HiC_analysis/extracting_loop_counts/correction_factors/{sample}.pk")

if __name__=="__main__":
    run_precompute()