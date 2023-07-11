# folder containing the preprocessing scripts for allelic imbalance calling for Hi-C data
These scripts are for reference only.

The process for calling allelic imbalance involves the following steps:

- Call genotypes from the Hi-C data. These are more complete than the genotype array data (genotype_calling.sh).
- Phase the genotypes using the Hi-C data through the integrated phasing pipeline (phasing_step1.sh).
- Split the bam files for each sample, this needs first the Hi-C pro bam files to be sorted (prepare_bams_per_sample.sh) and the SNPs file from phasing (prepare_snps.sh and py). Finally we use run_split.sh to retrieve the reads that align to each allele of the patient.