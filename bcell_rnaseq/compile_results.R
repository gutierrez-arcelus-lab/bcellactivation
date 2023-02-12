library(tidyverse)

gene_tx <- "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.txToGene.tsv" %>%
    read_tsv()

meta <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/metadata.tsv" |>
    read_tsv()

ids <- paste(meta$sample_id, meta$stim, sep = "_")  

id_df <- tibble(id = ids) |>
    separate(id, c("sample_id", "stim"), sep = "_", remove = FALSE)

results_df <- file.path("results", ids, "quant.sf") |>
    setNames(ids) |>
    map_df(read_tsv, .id = "id")

transcript_df <- results_df |>
    left_join(id_df, by = "id") |>
    left_join(gene_tx, by = c("Name" = "transcript_id")) |>
    select(sample_id, stim, gene_id, gene_name, tx_id = Name, 
	   counts = NumReads, tpm = TPM)

gene_df <- transcript_df |>
    group_by(sample_id, stim, gene_id, gene_name) |>
    summarise_at(vars(counts, tpm), sum) |>
    ungroup()

write_tsv(transcript_df, "./quantifications_transcripts.tsv")
write_tsv(gene_df, "./quantifications_genes.tsv")



#

ids <- paste(meta$sample_id, meta$stim, sep = "_")

mapped_reads <- sprintf("./mapping/%s_Log.final.out", ids) |>
    setNames(ids) |>
    map_df(~read_tsv(., col_names = FALSE) |>
	   filter(grepl("Uniquely mapped reads", X1)) |>
	   mutate(stat = str_extract(X1, "number|%"),
		  value = parse_number(X2)) |>
	   select(stat, value), .id = "id") |>
    separate(id, c("sample_id", "stim"), sep = "_")

mapped_reads |> filter(stat == "%") |> arrange(value)

mapped_reads |> 
    filter(stat == "number") |>
    select(sample_id, stim, reads = value) |>
    write_tsv("./star_n_uniq_mapped_reads.tsv")
    


