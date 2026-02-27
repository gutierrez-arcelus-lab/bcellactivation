library(tidyverse)
library(glue)

# Transcript annotations
gtf <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v41.primary_assembly.annotation.gtf.gz" |>
    rtracklayer::import(feature.type = "transcript")

gene_tx <- 
    gtf |>
    as.data.frame() |>
    as_tibble() |>
    select(transcript_id, gene_id, gene_name) |> 
    mutate(gene_id = str_remove(gene_id, "\\.\\d+"))

# Individual libraries
sample_ids <- 
    "./data/metadata_qced.tsv" |>
    read_tsv() |>
    pull(sample_id)

salmon_gene_df <- 
    glue("./salmon_quant/{sample_ids}/quant.sf") |>
    setNames(sample_ids) |>
    map_dfr(~read_tsv(.) |> 
	    left_join(gene_tx, join_by(Name == transcript_id)) |>
	    group_by(gene_id, gene_name) |>
	    summarize_at(vars(NumReads, TPM), sum) |> 
	    ungroup(),
	    .id = "sample_id") |>
    select(sample_id, gene_id, gene_name, read_count = NumReads, tpm = TPM)

write_tsv(salmon_gene_df, "./data/salmon_gene.tsv")


# Pooled replicates
pooled_sample_ids <- 
    "./data/metadata_qced_pooled.tsv" |>
    read_tsv(col_types = "c--") |>
    pull(sample_id)

pooled_meta_df <- 
    tibble(sample_id = pooled_sample_ids) |>
    separate(sample_id, c("stim", "timep", "donor_id"), sep = "\\.", remove = FALSE)

pooled_salmon_gene_df <- 
    glue("./salmon_quant_pooled/{pooled_sample_ids}/quant.sf") |>
    setNames(pooled_sample_ids) |>
    map_dfr(~read_tsv(.) |> 
	    left_join(gene_tx, join_by(Name == transcript_id)) |>
	    group_by(gene_id, gene_name) |>
	    summarize_at(vars(NumReads, TPM), sum) |> 
	    ungroup(),
	    .id = "sample_id") |>
    left_join(pooled_meta_df, join_by(sample_id)) |>
    select(donor_id, stim, timep, gene_id, gene_name, read_count = NumReads, tpm = TPM)

write_tsv(pooled_salmon_gene_df, "./data/salmon_pooledreps_gene.tsv")
