#!/bin/bash --login
#$ -pe smp.pe 4
#$ -j y
#$ -o scratch/logs
#$ -t 1-1364

###### 
# this script is for reference only
# it takes the bam file for each Hi-C aligned using bwa-mem
# bwa-0.7.17/bwa mem -SP5M -t 8 
# one file per patient

INDEX=$((SGE_TASK_ID-1))
sample_ID=$((${INDEX}/22))
CHR_N=$((${INDEX}-${sample_ID}*22+1))
echo $CHR_N
echo $sample_ID

samples=("NRHV014X" "NRHV073" "NRHV079" "NRHV086" "NRHV121" "NRHV151" "NRHV168" "NRHV171" "NRHV238" "NRHV290" "NRHV295" "NRHV321" "NRHV322" "NRHV324" "NRHV325" "NRHV326" "NRHV332" "PSA4917" "PSA4918" "PSA4920" "PSA4941" "PSA4942" "PSA4943" "PSA4944" "PSA4945" "PSA4946" "PSA4950" "PSA4951" "PSA4954" "PSA4957" "PSA4958" "PSA4959" "PSA4960" "PSA4961" "PSA4962" "PSA4963" "PSA4967" "PSA4968" "PSA4969" "PSA5006" "PSA5007" "PSA5008" "PSA5009" "PSA5010" "PSA5012" "PSA5013" "PSA5014" "PSA5015" "PSA5017" "PSA5018" "PSA5019" "PSA5020" "PSA5021" "PSA5022" "PSA5023" "PSA5024" "PSA5025" "PSA5026" "PSA5036" "PSA5037" "PSA5039" "PSA5040")
sample=${samples[$sample_ID]}
echo $sample


# step 2
# 3.2 of glimpse. BCFTOOLS mpileup
# modified to exclude the snps that overlap restriction cut sites
static_bins=/mnt/jw01-aruk-home01/projects/psa_functional_genomics/SSc_combined_analysis/GLIMPSE/static_bins

ref=/mnt/jw01-aruk-home01/projects/psa_functional_genomics/SSc_combined_analysis/GLIMPSE/dataset

source activate /mnt/jw01-aruk-home01/projects/psa_functional_genomics/SSc_combined_analysis/GLIMPSE_ENV

BAM=/net/scratch2/mdefscs4/HiC_variant_calling/${sample}/${sample}_merged_sorted.bam
REFGEN=/mnt/jw01-aruk-home01/projects/psa_functional_genomics/NEW_references/refseq_genome_indexes/GRCh38_no_alt_plus_hs38d1_analysis.fasta.gz
VCF=$ref/chr${CHR_N}.sites.vcf.gz
TSV=$ref/chr${CHR_N}.sites.tsv.gz
OUT=/net/scratch2/mdefscs4/HiC_variant_calling/${sample}/glimpse_chr${CHR_N}.vcf.gz
EXCL=/mnt/jw01-aruk-home01/projects/functional_genomics/bin/HICUP_0.7.4/arima_digest_refseq_noalt/cutsites_complement_arima_hg38_cleaned.bed.gz

bcftools mpileup -f ${REFGEN} -R ${EXCL} -I -E -a 'FORMAT/DP' -T ${VCF} -t chr${CHR_N} ${BAM} -Ou --ignore-RG --threads 4 | \
bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT} --threads 4

bcftools index -f ${OUT}

# step 3
# 5.1 of glimpse

VCF=/net/scratch2/mdefscs4/HiC_variant_calling/${sample}/glimpse_chr${CHR_N}.vcf.gz
REF=$ref/chr${CHR_N}.bcf
MAP=/mnt/jw01-aruk-home01/projects/psa_functional_genomics/SSc_combined_analysis/GLIMPSE/maps/genetic_maps.b38/chr${CHR_N}.b38.gmap.gz
while IFS="" read -r LINE || [ -n "$LINE" ];
do
 printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
 IRG=$(echo $LINE | cut -d" " -f3)
 ORG=$(echo $LINE | cut -d" " -f4)
 OUT=/net/scratch2/mdefscs4/HiC_variant_calling/${sample}/imputed_glimpse_chr${CHR_N}_${ID}.bcf
 GLIMPSE_phase --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
 bcftools index -f ${OUT}
done < $ref/chunks.chr${CHR_N}.txt

# step 4
# step 6.1 of glimpse
LST=/net/scratch2/mdefscs4/HiC_variant_calling/${sample}/list_chr${CHR_N}.txt
ls /net/scratch2/mdefscs4/HiC_variant_calling/${sample}/imputed_glimpse_chr${CHR_N}_*.bcf > ${LST}
OUT=/net/scratch2/mdefscs4/HiC_variant_calling/${sample}/merged_glimpse_chr${CHR_N}.bcf
GLIMPSE_ligate --input ${LST} --output $OUT
bcftools index -f ${OUT}
bcftools convert -O z -o /net/scratch2/mdefscs4/HiC_variant_calling/${sample}/merged_glimpse_chr${CHR_N}.vcf.gz ${OUT}
bcftools index -f /net/scratch2/mdefscs4/HiC_variant_calling/${sample}/merged_glimpse_chr${CHR_N}.vcf.gz

# step 5
# step 7.1 to get phased haplotypes

VCF=/net/scratch2/mdefscs4/HiC_variant_calling/${sample}/merged_glimpse_chr${CHR_N}.bcf
OUT=/net/scratch2/mdefscs4/HiC_variant_calling/${sample}/merged_glimpse_chr${CHR_N}_phased.bcf
GLIMPSE_sample --input ${VCF} --solve --output ${OUT}
bcftools index -f ${OUT}
bcftools convert -O z -o /net/scratch2/mdefscs4/HiC_variant_calling/${sample}/merged_glimpse_chr${CHR_N}_phased.vcf.gz ${OUT}
bcftools index -f /net/scratch2/mdefscs4/HiC_variant_calling/${sample}/merged_glimpse_chr${CHR_N}_phased.vcf.gz

# copy files that are needed

dest_folder=/mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/genotype_calling/called_genotypes
mkdir -p $dest_folder/$sample/per_chrom

cp /net/scratch2/mdefscs4/HiC_variant_calling/${sample}/merged_glimpse_chr${CHR_N}_phased.vcf.gz $dest_folder/$sample/per_chrom/
cp /net/scratch2/mdefscs4/HiC_variant_calling/${sample}/merged_glimpse_chr${CHR_N}_phased.bcf $dest_folder/$sample/per_chrom/
cp /net/scratch2/mdefscs4/HiC_variant_calling/${sample}/merged_glimpse_chr${CHR_N}_phased.vcf.gz.csi $dest_folder/$sample/per_chrom/
cp /net/scratch2/mdefscs4/HiC_variant_calling/${sample}/merged_glimpse_chr${CHR_N}_phased.bcf.csi $dest_folder/$sample/per_chrom/
source activate personal_software
tabix -p vcf $dest_folder/$sample/per_chrom/merged_glimpse_chr${CHR_N}_phased.vcf.gz


# rm files

rm /net/scratch2/mdefscs4/HiC_variant_calling/${sample}/*chr${CHR_N}_*.bcf*
rm /net/scratch2/mdefscs4/HiC_variant_calling/${sample}/*chr${CHR_N}.vcf*
rm /net/scratch2/mdefscs4/HiC_variant_calling/${sample}/*chr${CHR_N}_*.vcf*
rm /net/scratch2/mdefscs4/HiC_variant_calling/${sample}/*chr${CHR_N}.txt*
