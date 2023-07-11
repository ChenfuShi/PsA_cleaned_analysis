#!/bin/bash --login
#$ -pe smp.pe 2
#$ -j y
#$ -o scratch/logs
#$ -t 1-9

INDEX=$((SGE_TASK_ID-1))

cd /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_cleaned_analysis/QTL_analysis/HiC
module load apps/gcc/qtltools/1.3.1


numbers=(0 1 2 3 4 5 7 10 15)

for i in "${numbers[@]}"
do
    QTLtools cis --vcf /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/QTL_tools/genotypes_HiC_ATAC_merged_0.05.vcf.gz \
    --cov covariates/loop_covariates_CD8_${i}_${numbers[$INDEX]}.txt --permute 1000 \
    --normal \
    --bed loop_CD8_test.bed.gz --region chr1 --out output_permuted_test/loop_permuted_CD8_cov_${i}_${numbers[$INDEX]}.txt


    QTLtools cis --vcf /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/QTL_tools/genotypes_HiC_ATAC_merged_0.05.vcf.gz \
    --cov covariates/loop_covariates_CD4_${i}_${numbers[$INDEX]}.txt --permute 1000 \
    --normal \
    --bed loop_CD4_test.bed.gz --region chr1 --out output_permuted_test/loop_permuted_CD4_cov_${i}_${numbers[$INDEX]}.txt
done