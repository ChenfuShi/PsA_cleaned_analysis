import pandas as pd
import numpy as np
import pysam
import gzip
import os

downloaded_hic = pd.read_excel("/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_combined_analysis/metadata/downloaded_data.xlsx", sheet_name = "HiC")
downloaded_hic = downloaded_hic[downloaded_hic["condition"].isin(["patient", "synovium", "healthy"])]
downloaded_hic = downloaded_hic.drop_duplicates(subset = ["patient","cell_type"], keep = "first")

individuals = downloaded_hic["patient"].drop_duplicates().to_list()

chromosomes = list(range(1,23))

target_loc = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/allele_splitting/variant_files"

for sample in individuals:
    try:
        out_file = gzip.open(target_loc + f"/{sample}_SNPsplit.txt.gz","wt")
        out_file.write("ID\tChr\tPosition\tSNP value\tRef/SNP\n")
        for chrom in chromosomes:
            VCF_file = pysam.VariantFile(f'/mnt/iusers01/jw01/mdefscs4/psa_functional_genomics/PsA_combined_analysis/allele_specific_hic/genotype_calling/called_genotypes/{sample}/phased_per_chrom/HAPCUT_phased_chr{str(chrom)}.vcf.gz','r')
            res = VCF_file.fetch(contig = f"chr{str(chrom)}")
            for variant in res:
                s = variant.samples.items()[0][0]
                if variant.samples[s]["GT"] == (1,0) or variant.samples[s]["GT"] == (0,1):
                    out_file.write(f'{variant.id}\t{variant.chrom}\t{variant.pos}\t1\t{variant.alleles[variant.samples[s]["GT"][0]]}/{variant.alleles[variant.samples[s]["GT"][1]]}\n')
        out_file.close()
    except:
        print(sample, "failed, double check")


