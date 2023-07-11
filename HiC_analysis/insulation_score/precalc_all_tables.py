import pandas as pd
import numpy as np
from multiprocessing import Pool
import cooler
import cooltools
##############
# script to precalculate the insulation score table for samples


##############

import sys

base_dir = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_cleaned_analysis"

sample = sys.argv[1]

resolution = 25000
clr = cooler.Cooler(f'/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_HiC_analyzer/cool/{sample}.mcool::/resolutions/'+str(25000))

windows = [4*resolution, 5*resolution, 10*resolution]
insulation_table = cooltools.insulation(clr, windows, clr_weight_name=None,ignore_diags = 1,nproc=4)

insulation_table.to_csv(f"{base_dir}/HiC_analysis/insulation_score/ins_tables/{sample}.csv.gz")