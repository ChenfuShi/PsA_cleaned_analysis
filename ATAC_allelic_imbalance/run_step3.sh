#!/bin/bash --login
#$ -pe smp.pe 8
#$ -j y
#$ -o scratch/logs
#$ -t 1-7
##$ -l mem256

INDEX=$((SGE_TASK_ID-1))

cd /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_cleaned_analysis/ATAC_allelic_imbalance
source activate /mnt/jw01-aruk-home01/projects/functional_genomics/bin/new_basic_software

python /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_cleaned_analysis/ATAC_allelic_imbalance/3_combined_p_vals.py $INDEX