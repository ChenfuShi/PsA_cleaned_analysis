#!/bin/bash --login
#$ -pe smp.pe 4
#$ -j y
#$ -o scratch/logs
##$ -l v100           # A 1-GPU request (v100 is just a shorter name for nvidia_v100)


#$ -t 1-1 

source activate /mnt/jw01-aruk-home01/projects/functional_genomics/bin/new_basic_software

cd /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/allele_splitting

python prepare_snps.py