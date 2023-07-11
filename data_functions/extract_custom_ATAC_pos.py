import pandas as pd
import pyBigWig
from pandarallel import pandarallel

pandarallel.initialize()
base_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_cleaned_analysis"

metadata_ATAC = pd.read_csv(f"{base_dir}/metadata/cleaned_ATAC_metadata.csv", index_col = 0)
ATAC_bw_loc = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/coverages/{}_ATAC/{}_ATAC_coverage.bw"
size_factors = pd.read_csv(f"{base_dir}/ATAC_seq_analysis/diffbind_size_factors.csv")
size_factors.columns = "ID size_factor".split()


def get_region(chrom, start, end):
    def _get_normalized_bw(row):
        cur_bw_loc = ATAC_bw_loc.format(row["id"],row["id"])
        # Open the original .bw file
        bw = pyBigWig.open(cur_bw_loc)

        # Get the normalization factor
        normalization_factor = float(size_factors.loc[size_factors["ID"] == row["id"],"size_factor"])

        res = bw.stats(chrom,start, end, type="sum")

        res = res[0]/normalization_factor

        # Close the files
        bw.close()
        return pd.Series({"id":row["id"], "proper_name":row["proper_name"], "peak_height" : res})

    return metadata_ATAC.apply(_get_normalized_bw, axis = 1)