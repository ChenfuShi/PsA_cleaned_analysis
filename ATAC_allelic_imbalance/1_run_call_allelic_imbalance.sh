#!/bin/bash --login
#$ -pe smp.pe 2
#$ -j y
#$ -o scratch/logs
#$ -t 1-2772

# script that runs the bcftools to call allelic imbalance from bam files of ATAC-seq
# bam files are generated using the ATAC-seq pipelines

# note that we cannot provide raw genotypes for the patients, so this script is here only as a reference and it cannot be run

INDEX=$((SGE_TASK_ID-1))
sample_ID=$((${INDEX}/22))
CHR_N=$((${INDEX}-${sample_ID}*22+1))
echo $CHR_N
echo $sample_ID

samples=(NRHV151_CD4 NRHV151_CD8 NRHV321_CD4 PsA4920_CD4 PsA4945_CD8 PsA4948_CD8 PsA4953_CD8 PsA4960_CD8 NRHV321_CD8 PsA4944_CD4 PsA4947_CD8 PsA4949_CD4 PsA4954_CD8 PsA4963_CD8 NRHV238_CD4 nrhv238_CD8_ PsA4950_CD4 PsA4961_CD8 PsA4964_CD8 PsA5012_CD4 NRHV326_CD8 PsA4920_CD8 PsA4952_CD8 PsA4963_CD4 PsA5008_CD4 PsA5012_CD8 PsA4941CD8 PsA4944CD8 PSA4946CD8 PsA4952CD4 PsA4954CD4 PsA4962CD8 PSA4964CD4 PsA4965CD8 PSA4969CD8 PsA4969CD8_SF PSA5008CD8 PsA5014CD4 PsA5014CD8 PsA5019CD8 NRHV171CD4 NRHV171CD8 NRHV324CD4 PsA4950CD8 PsA4958CD4 PsA4958CD8 PsA5009CD4 PsA5010CD4 PsA5010CD8 PsA5018CD4 PsA5018CD8 PsA5023CD8 PsA5026CD4 PsA5026CD8 PsA4956CD4 PsA4956CD8 PsA4968CD8 PsA5007CD8 PsA5009CD8 PsA5013CD4 PsA5013CD8 PsA5015CD8 PsA5020CD8 PsA5021CD8 PsA5022CD8 PsA5024CD8 PsA5036CD4 PsA5036CD8 NRHV086CD4 NRHV086CD8 NRHV121CD4 NRHV121CD8 NRHV322CD4 NRHV322CD8 NRHV324CD8 PSA4947CD4 PsA4957CD4 PsA4959CD4 PsA5006CD8 PsA5019CD4 PsA5021CD4 NRHV073CD4 NRHV073CD8 NRHV079CD4 NRHV079CD8 NRHV168CD4 NRHV290CD4 NRHV295CD4 NRHV295CD8 NRHV332CD8 PSA5040CD8 PSA5040CD8SF NRHV168_CD8 NRHV290_CD8 NRHV326CD4 NRHV332_CD4 PSA4945_CD4 PSA4949_CD8 psa4951_CD4 psa4951_CD8 PSA4959_CD8 PSA5017_CD4 PSA5017_CD8 PSA5025CD8SF PSA5037_CD4 PSA5039_CD8 NRHV014XCD4 NRHV014XCD8 NRHV325_CD4 NRHV325_CD8 NRSF055CD4SF NRSF055CD8SF NRSF056CD8SF PSA4951CD4SF PSA4958CD8SF PSA4967CD4 PSA4968CD4 PSA5006CD4 PSA5007CD4 PSA4942CD8SF PSA4955CD8 PSA4958CD4SF PSA4962CD4 PSA4966CD8 PSA5037CD8 PSA5040CD4)
sample=${samples[$sample_ID]}
echo $sample

cd /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/ATAC_specific/called_bcf
mkdir -p $sample
cd $sample

BAM=/net/scratch2/mdefscs4/master_ATAC_ChIP_analyzer/clean_alignments_alelle/${sample}_ATAC/${sample}_ATAC_align_filtered_macs2.bam
REFGEN=/mnt/jw01-aruk-home01/projects/psa_functional_genomics/NEW_references/refseq_genome_indexes/GRCh38_no_alt_plus_hs38d1_analysis.fasta.gz
VCF=$ref/chr${CHR_N}.sites.vcf.gz
TSV=$ref/chr${CHR_N}.sites.tsv.gz
OUT=/net/scratch2/mdefscs4/master_ATAC_ChIP_analyzer/do_allele_call/${sample}
mkdir -p ${OUT}

bcftools mpileup -q 10 -f ${REFGEN} -I -E -a 'FORMAT/DP,FORMAT/AD' -T ${VCF} -t chr${CHR_N} ${BAM} -Ou --ignore-RG -d 10000 | \
 bcftools call -Aim -C alleles -T ${TSV} -Oz --threads 2 -o ${OUT}/chr${CHR_N}.called.vcf.gz

bcftools view -i 'MIN(FMT/DP)>10' ${OUT}/chr${CHR_N}.called.vcf.gz -Oz -o ${OUT}/chr${CHR_N}.called_filtered.vcf.gz


bcftools index ${OUT}/chr${CHR_N}.called_filtered.vcf.gz

bcftools annotate -a ~/psa_functional_genomics/NEW_references/variants/gatk_references/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
 -c ID -o /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/ATAC_specific/called_bcf/${sample}/${sample}_chr${CHR_N}_filt_annotated.vcf.gz \
 ${OUT}/chr${CHR_N}.called_filtered.vcf.gz

bcftools index /mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/ATAC_specific/called_bcf/${sample}/${sample}_chr${CHR_N}_filt_annotated.vcf.gz