# ==============================================================================
# Description:  Prepares ATAC-seq peaks and GWAS summary statistics for stratified 
#               LD Score Regression (S-LDSC). This involves creating configuration 
#               files for LDSC, using liftOver to convert GRCh38 peaks to GRCh37 
#               (required for standard LDSC), extracting condition-specific peak 
#               sets, and formatting the list of GWAS traits to be tested.
# Input:        1. ../1-processing/.../consensus_peaks.mLb.clN.featureCounts.txt
#               2. ../2-differential_peaks/results/*_24vsunst_24.tsv
#               3. /reference_databases/.../hg38ToHg19.over.chain.gz (LiftOver chain)
#               4. /lab-share/.../sumstats_Description_070121.xlsx (GWAS metadata)
# Output:       Peak coordinates, continuous annotations (.ldcts), and curated 
#               lists of GWAS summary statistics for downstream Slurm arrays.
# ==============================================================================

library(tidyverse)
library(glue)
library(readxl)

# ------------------------------------------------------------------------------
# 1. Environment & LDSC Configuration Setup
# ------------------------------------------------------------------------------
if (!file.exists("data")) dir.create("data")
if (!file.exists("data/ldscores")) dir.create("data/ldscores")
if (!file.exists("results")) dir.create("results")

stim_order <- c("IL4", "TLR7", "BCR", "DN2")

# Generate an array specification file for Slurm.
# LDSC annotates peaks chromosome by chromosome. This creates a grid mapping 
# each condition (and the background control) across chromosomes 1-22.
expand_grid(chr = 1:22, 
	    set = c(stim_order, "control")) |>
    mutate(set = fct_inorder(set)) |>
    arrange(chr, set) |>
    write_tsv("./data/array_spec.tsv", col_names = FALSE)

# Write the Continuous Trait Scores (CTS) file.
# This file tells the LDSC partitioned heritability algorithm which specific 
# annotation set to compare against the baseline control.
tibble(set = stim_order) |>
    mutate(annot = glue("data/ldscores/{set}."),
	   control = "data/ldscores/control.") |>
    unite("ldscores", c(annot, control), sep = ",") |>
    write_tsv("./data/peaksets.ldcts", col_names = FALSE)

# ------------------------------------------------------------------------------
# 2. Peak Import & Coordinate Liftover (GRCh38 -> GRCh37)
# ------------------------------------------------------------------------------
# Import the master consensus ATAC peaks (currently in GRCh38)
consensus_peaks <- 
    "../1-processing/results/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.featureCounts.txt" |>
    read_table(skip = 1) |>
    select(chrom = Chr, start = Start, end = End, peak_id = Geneid) |>
    filter(chrom %in% paste0("chr", c(1:22))) |>
    mutate(chrom = factor(chrom, levels = paste0("chr", c(1:22)))) |>
    arrange(chrom, start, end)

# Standard LDSC baseline models rely on GRCh37 (hg19) coordinates. 
# We use UCSC liftOver to map our GRCh38 peaks to GRCh37.
bed38_file <- "./data/hg38.bed"
bed19_file <- "./data/hg19.bed"
fail_file <- "./data/failTolift.txt"
# NOTE FOR EXTERNAL USERS: provide the appropriate path to the UCSC chain file 
chain_file <- "/reference_databases/ReferenceGenome/liftover_chain/hg38/hg38ToHg19.over.chain.gz"

write_tsv(consensus_peaks, bed38_file, col_names = FALSE)

# Execute system command to run liftOver
glue("liftOver {bed38_file} {chain_file} {bed19_file} {fail_file}") |>
    system()

# Read the successfully converted GRCh37 peaks
bedlift <- read_tsv(bed19_file, col_names = c("chrom", "start", "end", "peak_id"))

# Format the coordinates into a 'gene-coord-file' format required by make_annot.py
bedlift |>
    select(GENE = peak_id, CHR = chrom, START = start, END = end) |>
    write_tsv("./data/peak_coord.txt")

# Save all successfully converted peak IDs to serve as the background 'control'
write_lines(bedlift$peak_id, "./data/control.txt")

# ------------------------------------------------------------------------------
# 3. Differentially Accessible Peak Extraction
# ------------------------------------------------------------------------------
# Parse the DESeq2 outputs to identify peaks that become significantly more 
# accessible (padj <= 0.01, log2FC >= 1) upon stimulation.
da_files <- list.files("../2-differential_peaks/results", pattern = ".+vs.+\\.tsv$", full.names = TRUE)
da_files <- setNames(da_files, basename(da_files) |> str_remove("\\.tsv"))
da_files <- keep(da_files, grepl("vsunst_24", names(da_files)))

da_data <-
    map_dfr(da_files, read_tsv, .id = "comparison") |>
    filter(!is.na(padj)) |>
    separate(comparison, c("stim1", "stim2"), sep = "vs") |>
    select(stim = stim1, peak_id = interval, chrom = Chr, log2fc = log2FoldChange, padj) |>
    mutate(stim = str_remove(stim, "_24")) |>
    # Inner join with the lifted-over peaks to ensure we only test peaks 
    # that successfully mapped to GRCh37
    inner_join(bedlift, join_by(chrom, peak_id))

# Split the significant peak IDs by condition and export them to text files
set_list <- 
    da_data |> 
    filter(padj <= 0.01, log2fc >= 1) |>
    arrange(chrom, start, end) |>
    select(stim, peak_id) |>
    {function(x) split(x, x$stim)}() |>
    map(~pull(., peak_id))

walk2(set_list, names(set_list), 
     ~write_lines(.x, glue("./data/{.y}.txt")))

# ------------------------------------------------------------------------------
# 4. GWAS Trait Curation & Formatting
# ------------------------------------------------------------------------------
# Parse metadata to link disease traits (e.g., "Asthma") to their specific LDSC 
# summary statistic file paths.
#
# The summary statistics and reference files were obtained from:
# https://zenodo.org/records/7768714

# A) Original institutional LDSC summary stats
## Original summary stats
sumstats_list <- 
    "./data/ldsc_data/sumstats_Description_070121.xlsx" |>
    read_excel() |>
    janitor::clean_names() |>
    select(trait_name, trait_identifier, reference)

sumstats_files <- 
    list.files("./data/ldsc_data/v4/sumstats/")

# Define the specific immune, autoimmune, and control traits of interest for this paper
traits <- 
    c("Height", "HDL", "LDL", "Crohn's Disease", "Rheumatoid Arthritis", "Schizophrenia",
      "Type 1 Diabetes", "Type 2 Diabetes", "Ulcerative Colitis", "IBD", "Multiple Sclerosis",
      "Celiac Disease", "Primary Biliary Cirrhosis", "Systemic Lupus Erythematosus",
      "Eczema", "Asthma", "Adult Onset Asthma", "Child Onset Asthma", "Psoriasis", 
      "MDD", "Cancer", "SARS-CoV-2 infection", "Covid 19 Vaccination")

sumstats_df <- 
    sumstats_list |>
    mutate(trait_name = str_remove(trait_name, " \\(.+\\)$")) |>
    mutate(trait_name = case_when(tolower(trait_name) == "major depressive disorder" ~ "MDD",
				  trait_name == "High Density Lipoprotein" ~ "HDL",
				  TRUE ~ trait_name)) |>
    filter(trait_name %in% traits,
	   trait_identifier %in% str_remove(sumstats_files, "\\.sumstats\\.gz$")) |>
    arrange(trait_name)

# B) Supplementary 107 independent traits database
sumstats_107_list <-
    "./data/ldsc_data/v4/sumstats_107/traits_indep107.xlsx" |>
    read_excel() |>
    filter(is.na(`...6`)) |>
    janitor::clean_names() |>
    select(trait, file_name, paper_link) |> 
    mutate(trait = str_remove(trait, "T1D - "),
	   trait = str_replace(trait, "Lupus\\s+\\(SLE\\)", "Systemic Lupus Erythematosus"),
	   trait = str_replace(trait, "Celiac disease", "Celiac Disease"),
	   trait = str_remove(trait, " \\(.+\\)$")) |>
    select(trait_name = trait, trait_identifier = file_name, reference = paper_link)

sumstats_107_files <- 
    "./data/ldsc_data/v4/sumstats_107/" |>
    list.files(pattern = "\\.sumstats\\.gz$")

sumstats_107_df <- 
    sumstats_107_list |>
    filter(trait_name %in% traits,
	   trait_identifier %in% str_remove(sumstats_107_files, "\\.sumstats\\.gz$"))

# C) Merge and resolve duplicate traits between the two databases
# We explicitly drop duplicate representations of specific diseases from the 107 list 
# to prioritize the standard institutional summaries.
sumstats_all <- 
    bind_rows(sumstats = sumstats_df, sumstats_107 = sumstats_107_df, .id = "directory") |>
    mutate(trait_name = case_when(trait_name == "Adult Onset Asthma" ~ "Asthma (Adult Onset)",
				  trait_name == "Child Onset Asthma" ~ "Asthma (Child Onset)",
				  TRUE ~ trait_name)) |>
    arrange(trait_name, directory) |>
    mutate(trait_name = case_when(grepl("Covid19_Infection", trait_identifier) ~ "Covid19 Infection",
				  grepl("Covid19_Vaccination", trait_identifier) ~ "Covid19 Vaccination",
				  TRUE ~ trait_name)) |>
    filter(trait_name != "Covid19 Vaccination",
	   !(directory == "sumstats_107" & trait_name == "Crohn's Disease"),
	   !(directory == "sumstats_107" & trait_name == "Celiac Disease"),
	   !(directory == "sumstats_107" & trait_name == "HDL"),
	   !(directory == "sumstats_107" & trait_name == "LDL"),
	   !(directory == "sumstats_107" & trait_name == "IBD"),
	   !(directory == "sumstats_107" & trait_name == "Primary Biliary Cirrhosis"),
	   !(directory == "sumstats_107" & trait_name == "Systemic Lupus Erythematosus"))

# Export the final, clean list of trait paths for the downstream LDSC array job
write_tsv(sumstats_all, "./data/traits.txt", col_names = FALSE)
