# ==============================================================================
# Description:  Imports transcript-level quantifications from Salmon, maps 
#               transcripts to their parent genes using GENCODE annotations, 
#               and aggregates the data into gene-level read counts and TPMs.
# Input:        1. gencode.v38.primary_assembly.annotation.gtf.gz (Annotations)
#               2. metadata.tsv (Sample IDs)
#               3. quant.sf files (Salmon output per sample)
# Output:       1. expression_pooled_reps_transcriptlevel.tsv
#               2. expression_pooled_reps.tsv (Gene-level aggregated data)
# ==============================================================================
library(tidyverse)
library(tximport)

# ------------------------------------------------------------------------------
# 1. Transcript-to-Gene Mapping
# ------------------------------------------------------------------------------
# Import the transcript features from the GENCODE GTF
gtf <- 
    "./data/gencode.v38.primary_assembly.annotation.gtf.gz" |>
    rtracklayer::import(feature.type = "transcript")

# Create a mapping table connecting transcript IDs to gene IDs and names
gene_tx <- 
    gtf |>
    as_tibble() |>
    select(transcript_id, gene_id, gene_name) |> 
    mutate(gene_id = str_remove(gene_id, "\\.\\d+"))

# ------------------------------------------------------------------------------
# 2. Import Salmon Quantifications
# ------------------------------------------------------------------------------
# Read the sample metadata (expecting headerless: sample_id, fastq_paths)
meta <- 
    read_tsv("./data/metadata.tsv", col_names = c("sample_id", "fq1", "fq2"))

salmon_files <- 
    sprintf("./results/%s/quant.sf", meta$sample_id) |>
    setNames(meta$sample_id)

# Iterate over all samples, load their quant.sf files, and merge with annotations
salmon_df <- 
    salmon_files |>
    map_dfr(~read_tsv(.) |> 
	    left_join(gene_tx, join_by(Name == transcript_id)),
	    .id = "sample_id") |>
    select(sample_id, gene_id, gene_name, tx_id = Name, read_count = NumReads, tpm = TPM)

# ------------------------------------------------------------------------------
# 3. Aggregate to Gene Level
# ------------------------------------------------------------------------------
# Sum the transcript-level reads and TPMs to calculate total gene-level expression
salmon_gene <- 
    salmon_df |>
    group_by(sample_id, gene_id, gene_name) |>
    summarise_at(vars(read_count, tpm), sum) |> 
    ungroup()

# ------------------------------------------------------------------------------
# 4. Export Compiled Tables
# ------------------------------------------------------------------------------
write_tsv(salmon_df, "./results/expression_pooled_reps_transcriptlevel.tsv")
write_tsv(salmon_gene, "./results/expression_pooled_reps.tsv")

# ------------------------------------------------------------------------------
# 5. Import with tximport
# ------------------------------------------------------------------------------
txi <- 
    salmon_files |>
    tximport(type = "salmon", tx2gene = select(gene_tx, transcript_id, gene_id))

# Save data
write_rds(txi, "./results/gene_expression_txi.rds")
