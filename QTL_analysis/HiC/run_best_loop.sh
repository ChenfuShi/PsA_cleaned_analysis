#!/bin/bash --login
#$ -pe smp.pe 2
#$ -j y
#$ -o scratch/logs
#$ -t 1-23

INDEX=$((SGE_TASK_ID-1))

cd /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_cleaned_analysis/QTL_analysis/HiC
module load apps/gcc/qtltools/1.3.1


QTLtools cis --vcf /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/QTL_tools/genotypes_HiC_ATAC_merged_0.05.vcf.gz \
--cov covariates/loop_covariates_CD8_0_0.txt --nominal 0.01 \
--normal --chunk $INDEX 22 \
--bed loop_CD8_test.bed.gz --out output_final/loop_nominal_CD8_0_0_c$INDEX.txt


QTLtools cis --vcf /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/QTL_tools/genotypes_HiC_ATAC_merged_0.05.vcf.gz \
--cov covariates/loop_covariates_CD4_0_0.txt --nominal 0.01 \
--normal --chunk $INDEX 22 \
--bed loop_CD4_test.bed.gz --out output_final/loop_nominal_CD4_0_0_c$INDEX.txt


QTLtools cis --vcf /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/QTL_tools/genotypes_HiC_ATAC_merged_0.05.vcf.gz \
--cov covariates/loop_covariates_CD8_0_0.txt --permute 1000 \
--normal --chunk $INDEX 22 \
--bed loop_CD8_test.bed.gz --out output_final/loop_permuted_CD8_0_0_c$INDEX.txt


QTLtools cis --vcf /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/QTL_tools/genotypes_HiC_ATAC_merged_0.05.vcf.gz \
--cov covariates/loop_covariates_CD4_0_0.txt --permute 1000 \
--normal --chunk $INDEX 22 \
--bed loop_CD4_test.bed.gz --out output_final/loop_permuted_CD4_0_0_c$INDEX.txt