# Multi-omics analysis in primary T cells elucidates mechanisms behind disease associated genetic loci

## Introduction
This repository contains the scripts to reproduce the analysis presented in the manuscript.
The scripts will often use paths as set up to work in the environment of the university of manchester cluster. These will need to be changed when you download the preprocessed data and set up the scripts in your environment.

Please find all pre-processed data in http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/.
Note however that bam files will not be available as the reads information is considered protected access due to being able to call genotypes from the reads.

The repository on github will not include a lot of data due to size limits by github. You can find the full folder with all the files on http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/PsA_cleaned_analysis.

Two notebooks that might be useful:  
print all QTL and allelic imbalance information for a particular SNP: [access on colab](https://colab.research.google.com/github/ChenfuShi/PsA_cleaned_analysis/blob/main/integration_analysis/Everything_printer_public.ipynb)  
Identify the closes SNPs for which we have pregenerated a Hi-C correlation map: [access on colab](https://colab.research.google.com/github/ChenfuShi/PsA_cleaned_analysis/blob/main/integration_analysis/find_SNP.ipynb)  

As Genotypes and raw reads for patients are considered personal information these cannot be released. If you require these please get in contact with gisela.orozco[at]manchester.ac.uk

## Citation

If you use our data or use our code to run your analysis please cite the manuscript:

**Multi-omics analysis in primary T cells elucidates mechanisms behind disease associated genetic loci**  
Chenfu Shi, Danyun Zhao, Stefano Rossi, Antonios Frantzeskos, James Ding, Carlo Ferrazzano, Charlotte Wynn, Ryan Hum, Ellie Richards, Muskan Gupta, Chuan Fu Yap, Darren Plant, Richard Grencis, Paul Martin, Antony Adamson, Stephen Eyre, John Bowes, Anne Barton, Pauline Ho, Magnus Rattray, Gisela Orozco  
medRxiv 2023.07.19.23292550; doi: https://doi.org/10.1101/2023.07.19.23292550

## Requirements
The analysis includes a mix of Python and R scripts.
Python is run through a conda environment which contains most of the other tools necessary the analysis. The full environment definition is available in conda_env.yml.
In the scripts this is activated as new_basic_software, so if you want to use that script you need to change the conda activation to the name of the environment you have installed.


## Data Preprocessing
Data has been preprocessed according to the manuscript's methods. Preprocessing for ATAC-seq and Hi-C data is carried out using the pipeline available at https://github.com/ChenfuShi/hic_master_pipeline and https://github.com/ChenfuShi/ATAC_ChIP_pipeline.

## Useful Files
You can find the following files to be particularly useful if you want to do further analysis:  
QTL tables:

[eQTL nominal CD4+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/RNA/RNA_nominal_CD4_merged.txt)  
[eQTL nominal CD8+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/RNA/RNA_nominal_CD8_merged.txt)  
[eQTL permutation CD4+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/RNA/RNA_permuted_CD4_FDR.txt)  
[eQTL permutation CD8+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/RNA/RNA_permuted_CD8_FDR.txt)  


[caQTL nominal CD4+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/ATAC/ATAC_nominal_CD4_merged.txt)  
[caQTL nominal CD8+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/ATAC/ATAC_nominal_CD8_merged.txt)  
[caQTL permutation CD4+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/ATAC/ATAC_permuted_CD4_FDR.txt)  
[caQTL permutation CD8+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/ATAC/ATAC_permuted_CD8_FDR.txt)  


[insQTL nominal CD4+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/HiC/ins_nominal_CD4_merged.txt)  
[insQTL nominal CD8+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/HiC/ins_nominal_CD8_merged.txt)  
[insQTL permutation CD4+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/HiC/ins_permuted_CD4_FDR.txt)  
[insQTL permutation CD8+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/HiC/ins_permuted_CD8_FDR.txt)  


[loopQTL nominal CD4+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/HiC/loop_nominal_CD4_merged.txt)  
[loopQTL nominal CD8+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/HiC/loop_nominal_CD8_merged.txt)  
[loopQTL permutation CD4+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/HiC/loop_permuted_CD4_FDR.txt)  
[loopQTL permutation CD8+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/QTLS/HiC/loop_permuted_CD8_FDR.txt)  


Allelic imbalance tables:

[loop allelic imbalance CD4+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/hic/allelic_imbalance/allelic_imbalance_CD4_apeglm_results.csv)  
[loop allelic imbalance CD8+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/hic/allelic_imbalance/allelic_imbalance_CD8_apeglm_results.csv)  
[loop allelic imbalance merged](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/hic/allelic_imbalance/allelic_imbalance_ALL_apeglm_results.csv)  


[ATAC allelic imbalance CD4+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/atac/allelic_imbalance/ATAC_CD4_allelic_imbalance_with_betabinom.csv.gz)  
[ATAC allelic imbalance CD8+](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/atac/allelic_imbalance/ATAC_CD8_allelic_imbalance_with_betabinom.csv.gz)  
[ATAC allelic imbalance merged](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/atac/allelic_imbalance/ATAC_ALL_allelic_imbalance_with_betabinom.csv.gz)  


pregenerated correlation Hi-C maps for all genes, chromatin accessibility and highly significant loop and insulation QTLs:  
[webpage](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/PsA_output_hic_plots/main.html)

[samples metadata tables](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/metadata/)

[RNA-seq counts for all samples](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/rna/RNA_normalized_counts.csv)  
[ATAC-seq counts for all samples](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/atac/merged/ATAC_DESeq2_quantile_normalized_counts.csv)  


[Merged Hi-C juicebox files for visualization at up to 1kb resolution](http://bartzabel.ls.manchester.ac.uk/orozcolab/SNP2Mechanism/hic/merged/)

## Directory Listing

- ATAC_allelic_imbalance: scripts to call chromatin accessibility allelic imbalance, annotate the resulting tables and annotate GWAS studies with the results. In this folder you can also find the precomputed allelic imbalance (ATAC_allelic_imbalance/combined_p_vals_files/all_SNPs_all.csv) for all the SNPs that overlap chromatin accessibility regions in T cells. 
- ATAC_seq_analysis: Diffbind preprocessing of the data.
- data_functions: two helper functions for other scripts.
- HiC_allelic_imbalance: scripts to call loop allelic imbalance, including genotype calling and phasing from reads. In this folder you can also find the precomputed allelic imbalance loops (HiC_allelic_imbalance/output_dataframe_CD8.csv and HiC_allelic_imbalance/output_dataframe_CD4.csv)
- HiC_analysis: scripts to call mustache loops, extract loop counts, differential loop calling, insulation score, differential insulation score, outlier analysis from cell lines.
- integration_analysis: scripts that combine different omics together. Visualizations for changes in chromatin conformation with genotype, ATAC-peaks and gene expression. Correlation between loop strength with chromatin accessibility and gene expression. Correlation between insulation score and gene expression. scripts that allow the plotting of the data for specific regions (integration_analysis/plotter.ipynb) and printing of all the results from a specific set of SNPs (integration_analysis/everything_printer.ipynb).
- metadata: metadata files and GWAS snps with LD.
- QTL_analysis: All QTL analysis for the different omics, including files to annotate GWAS studies and compute overlaps between different QTLs and allelic imbalance.
- RNA_seq_analysis: DESeq2 preprossing of the data.
