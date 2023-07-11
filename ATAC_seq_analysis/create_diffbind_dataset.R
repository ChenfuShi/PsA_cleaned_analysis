#!/usr/bin/env Rscript

################
# This R script pregenerates the diffbind dataset across the whole dataset
# This is necessary because this step is very slow as it has to access the bam files for all the samples

# Diffbind is run with generally standard settings
# peak size has been set to 500bp as this improved results significantly
# a minimum of 5 samples are needed for a peak to be considered in the consensus set

################

library("tidyverse")
library("DiffBind")
library("readxl")
library("pheatmap")
options(bitmapType="cairo")

dataset_info_file = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/PsA_cleaned_analysis/metadata/cleaned_ATAC_metadata.csv"
dataset_info = read.csv(dataset_info_file)
dataset_peaks_location = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/macs2"
dataset_alignment_location = "/mnt/jw01-aruk-home01/projects/psa_functional_genomics/master_ATAC_ChIP_analyzer/clean_alignments"

dataset_info <- dataset_info %>% filter(condition %in% c("healthy", "patient", "synovium"))

data_for_diffbind <- dataset_info %>% select(id, patient, cell_type, condition) %>% mutate(folder = paste0(id, "_ATAC"))
data_for_diffbind <- data_for_diffbind %>% mutate(Peaks = paste0(dataset_peaks_location, "/", folder, "/", folder, "_peaks_nosex.narrowPeak"), 
bamReads = paste0(dataset_alignment_location, "/", folder, "/", folder, "_align_filtered_macs2.bam"), PeakCaller = "narrow")
data_for_diffbind = rename(data_for_diffbind, "id" = "SampleID", "condition" = "Condition", "cell_type" = "Tissue")

data_object <- dba(sampleSheet=data_for_diffbind, minOverlap = 5)
data_object <- dba.blacklist(data_object, blacklist=DBA_BLACKLIST_HG38, cores=16)

data_object$config$cores = 16
data_object <- dba.count(data_object, summits = 250)

save(data_object, file = "diffbind_object.Rdata")