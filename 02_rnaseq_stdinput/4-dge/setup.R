# ==============================================================================
# Description:  Prepares metadata for differential gene expression (DGE) analysis.
#               To allow direct comparison of fold-changes between standard-input 
#               and low-input RNA-seq (per reviewer request), samples must be 
#               re-quantified using the GENCODE v38 annotation. This script 
#               formats the metadata and pools technical replicates per donor 
#               into comma-separated lists for dynamic Salmon quantification.
# Input:        ../1-mapping/data/metadata.tsv
# Output:       ./data/metadata.tsv (Pooled sample sheet for Salmon)
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# 1. Environment Setup
# ------------------------------------------------------------------------------

# Ensure required subdirectories exist
if (!file.exists("data")) dir.create("data")
if (!file.exists("results")) dir.create("results")

# ------------------------------------------------------------------------------
# 2. Metadata Import & Reshaping
# ------------------------------------------------------------------------------

# Read the original mapping metadata
meta <- read_tsv("../1-mapping/data/metadata.tsv")

# Process the metadata to pool technical replicates.
meta_pooled <- 
    meta |>
    separate(sample_id, c("stim", "donor_id", "replic_id"), sep = "_") |>
    mutate(stim = recode(stim,
			 "unstday0" = "unst.0", 
			 "BCR" = "BCR.24", 
			 "TLR7" = "TLR7.24",
			 "DN2" = "DN2.24"),
	   stim = factor(stim, levels = c("unst.0", "TLR7.24", "BCR.24", "DN2.24"))) |>
    arrange(donor_id, replic_id, stim) |>
    group_by(donor_id, stim) |>
    summarise_at(vars(R1:R2), ~paste(., collapse = ",")) |>
    ungroup() |>
    unite("sample_id", c(stim, donor_id), sep = ".") |>
    select(sample_id, R1, R2)
    
# ------------------------------------------------------------------------------
# 3. Export
# ------------------------------------------------------------------------------

# Write the newly formatted, pooled metadata sheet for the Slurm array
write_tsv(meta_pooled, "./data/metadata.tsv")
