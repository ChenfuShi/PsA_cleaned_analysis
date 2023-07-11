#!/bin/bash --login
#$ -pe smp.pe 4
#$ -j y
#$ -o scratch/logs
##$ -l mem512
#$ -t 400-418 


###### 
# this script is for reference only
# this script takes the genotypes and the bwa-mem mapped file and uses the integrated phasing pipeline to phase the genotypes


INDEX=$((SGE_TASK_ID-1))
sample_ID=$((${INDEX}/22))
CHR_N=$((${INDEX}-${sample_ID}*22+1))
echo $CHR_N
echo $sample_ID

samples=("NRHV014X" "NRHV073" "NRHV079" "NRHV086" "NRHV121" "NRHV151" "NRHV168" "NRHV171" "NRHV238" "NRHV290" "NRHV295" "NRHV321" "NRHV322" "NRHV324" "NRHV325" "NRHV326" "NRHV332" "PSA4917" "PSA4918" "PSA4920" "PSA4941" "PSA4942" "PSA4943" "PSA4944" "PSA4945" "PSA4946" "PSA4950" "PSA4951" "PSA4954" "PSA4957" "PSA4958" "PSA4959" "PSA4960" "PSA4961" "PSA4962" "PSA4963" "PSA4967" "PSA4968" "PSA4969" "PSA5006" "PSA5007" "PSA5008" "PSA5009" "PSA5010" "PSA5012" "PSA5013" "PSA5014" "PSA5015" "PSA5017" "PSA5018" "PSA5019" "PSA5020" "PSA5021" "PSA5022" "PSA5023" "PSA5024" "PSA5025" "PSA5026" "PSA5036" "PSA5037" "PSA5039" "PSA5040")
sample=${samples[$sample_ID]}
echo $sample


# step 1
source activate /mnt/jw01-aruk-home01/projects/ra_challenge/Chris_Leung/Software/conda_hapcut2_1.3.2

mkdir -p /net/scratch2/mdefscs4/HiC_variant_calling/$sample
cd /net/scratch2/mdefscs4/HiC_variant_calling/$sample
mkdir -p chr$CHR_N
cd chr$CHR_N


cp /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/genotype_calling/called_genotypes/$sample/per_chrom/merged_glimpse_chr${CHR_N}_phased.vcf.gz ./
gunzip merged_glimpse_chr${CHR_N}_phased.vcf.gz


samtools quickcheck -v /net/scratch2/mdefscs4/HiC_variant_calling/temp_trim/${sample}/merged_sorted.bam && echo 'all ok' || exit

extractHAIRS --hic 1 --bam /net/scratch2/mdefscs4/HiC_variant_calling/temp_trim/${sample}/merged_sorted.bam \
--VCF merged_glimpse_chr${CHR_N}_phased.vcf \
--region chr$CHR_N \
--out /net/scratch2/mdefscs4/HiC_variant_calling/$sample/chr$CHR_N/extractHAIRS_fragment_file



# step 2
/mnt/jw01-aruk-home01/projects/ra_challenge/Chris_Leung/Software/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit \
-check --input-vcf merged_glimpse_chr${CHR_N}_phased.vcf \
-R /mnt/jw01-aruk-home01/projects/psa_functional_genomics/SSc_combined_analysis/shapeit2_ref/hg38_1k_genomes_chr${CHR_N}.hap.gz \
/mnt/jw01-aruk-home01/projects/psa_functional_genomics/SSc_combined_analysis/shapeit2_ref/hg38_1k_genomes_chr${CHR_N}.legend.gz \
/mnt/jw01-aruk-home01/projects/psa_functional_genomics/SSc_combined_analysis/shapeit2_ref/1k_genomes.sample \
--output-log shapeit_chr${CHR_N}.log

/mnt/jw01-aruk-home01/projects/ra_challenge/Chris_Leung/Software/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit \
--input-vcf merged_glimpse_chr${CHR_N}_phased.vcf \
-R /mnt/jw01-aruk-home01/projects/psa_functional_genomics/SSc_combined_analysis/shapeit2_ref/hg38_1k_genomes_chr${CHR_N}.hap.gz \
/mnt/jw01-aruk-home01/projects/psa_functional_genomics/SSc_combined_analysis/shapeit2_ref/hg38_1k_genomes_chr${CHR_N}.legend.gz \
/mnt/jw01-aruk-home01/projects/psa_functional_genomics/SSc_combined_analysis/shapeit2_ref/1k_genomes.sample \
-M /mnt/jw01-aruk-home01/projects/psa_functional_genomics/SSc_combined_analysis/shapeit2_ref/genetic_maps_hg38/genetic_map_chr${CHR_N}.txt \
--thread 4 \
--output-graph merged_shapeit_chr${CHR_N}_phased.graph 








# step 3
python2.7 /mnt/jw01-aruk-home01/projects/ra_challenge/Chris_Leung/Software/IntegratedPhasing/samplehaps.py \
merged_glimpse_chr${CHR_N} merged_shapeit_chr${CHR_N}_phased.graph 1000

python2.7 /mnt/jw01-aruk-home01/projects/ra_challenge/Chris_Leung/Software/IntegratedPhasing/encodereads.py \
merged_glimpse_chr${CHR_N}.hapsamples > step4_chr${CHR_N}.pseudo_reads






# step4
cat extractHAIRS_fragment_file step4_chr${CHR_N}.pseudo_reads > concat_chr${CHR_N}.fragments



HAPCUT2 --fragments concat_chr${CHR_N}.fragments  \
--VCF merged_glimpse_chr${CHR_N}_phased.vcf \
--output HAPCUT_integrated_phased \
--hic 1




# copying

target_dir=/mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/genotype_calling/called_genotypes/$sample

mkdir -p $target_dir/phased_per_chrom

conda deactivate
source activate personal_software


mv HAPCUT_integrated_phased.phased.VCF HAPCUT_phased_chr${CHR_N}.vcf
bgzip HAPCUT_phased_chr${CHR_N}.vcf
tabix -p vcf HAPCUT_phased_chr${CHR_N}.vcf.gz
cp HAPCUT_phased_chr${CHR_N}.vcf.gz* $target_dir/phased_per_chrom/

gzip HAPCUT_integrated_phased
mv HAPCUT_integrated_phased.gz HAPCUT_blocks_info_chr${CHR_N}.gz
cp HAPCUT_blocks_info_chr${CHR_N}.gz $target_dir/phased_per_chrom/

# rm all
cd
rm -r /net/scratch2/mdefscs4/HiC_variant_calling/$sample/chr$CHR_N
