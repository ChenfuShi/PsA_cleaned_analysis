#!/bin/bash --login
#$ -j y
#$ -o scratch/logs
##$ -l v100           # A 1-GPU request (v100 is just a shorter name for nvidia_v100)
#$ -pe smp.pe 3 # 8 CPU cores available to the host code
#$ -t 1-3

# script to run the python script to extract counts from the hic files for the merged files, as these have different names patterns

hic_loc=/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_HiC_analyzer/juiceboxes/
loops_file_CD4_8=/mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_cleaned_analysis/HiC_analysis/mustache_loops/merged_CD4_CD8_loops.bed
loops_file_CD8_8SF=/mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_cleaned_analysis/HiC_analysis/mustache_loops/merged_CD8SF_CD8_loops.bed
out_loc=/mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_cleaned_analysis/HiC_analysis/extracting_loop_counts/counts

INDEX=$((SGE_TASK_ID-1))
FILES=(big_merger_all_CD4_oct22
big_merger_all_CD8_oct22
big_merger_all_CD8_SF_oct22)

source activate /mnt/iusers01/jw01/mdefscs4/communal_software/new_basic_software

python /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_cleaned_analysis/HiC_analysis/extracting_loop_counts/retrieve_loop_counts.py \
-f $hic_loc/${FILES[$INDEX]}.hic \
-i $loops_file_CD4_8 \
-r 5000 \
-o $out_loc/merged_files/${FILES[$INDEX]}_counts_CD8_CD4


python /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_cleaned_analysis/HiC_analysis/extracting_loop_counts/retrieve_loop_counts.py \
-f $hic_loc/${FILES[$INDEX]}.hic \
-i $loops_file_CD8_8SF \
-r 5000 \
-o $out_loc/merged_files/${FILES[$INDEX]}_counts_CD8_CD8SF