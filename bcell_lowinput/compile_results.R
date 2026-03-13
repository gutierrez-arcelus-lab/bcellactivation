library(tidyverse)
library(glue)

# Transcript annotations
gtf <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf.gz" |>
    rtracklayer::import(feature.type = "transcript")

gene_tx <- 
    gtf |>
    as.data.frame() |>
    as_tibble() |>
    select(transcript_id, gene_id, gene_name) |> 
    mutate(gene_id = str_remove(gene_id, "\\.\\d+"))

# QC
meta_long <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_lowinput/metadata/description_hsapiens_qc.tsv" |>
    read_tsv() |>
    pivot_longer(fq1:fq2, names_to = "dummy", values_to = "fastq") |>
    separate_rows(fastq, sep = ",") |>
    mutate(fastq = basename(fastq),
           fastq = sub("\\.fq\\.gz", "", fastq),
	   lane = sub("^(\\d)_.+$", "\\1", fastq)) |>
    select(-barcode_seq, -dummy)

meta_df <- 
    meta_long |>
    filter(sample_id != "BLANK") |>
    distinct(plate, well, sample_id, stim, time) |>
    unite("sample_label", c(plate, well), sep = "_", remove = FALSE)
   
salmon_files <- 
    glue("./results/salmon/{meta_df$sample_label}/quant.sf") |>
    setNames(meta_df$sample_label)

salmon_gene_df <-
    salmon_files |>
    map_dfr(~read_tsv(.) |> 
	    left_join(gene_tx, join_by(Name == transcript_id)) |>
	    group_by(gene_id, gene_name) |>
	    summarize_at(vars(NumReads, TPM), sum) |> 
	    ungroup(),
	    .id = "sample_label") |>
    left_join(meta_df, join_by(sample_label)) |>
    select(sample_id, stim, time, gene_id, gene_name, read_count = NumReads, tpm = TPM)

write_tsv(salmon_gene_df, "./data/salmon_gene_quants.tsv")


#fastqc <-
#    "/temp_work/ch229163/fastq/lowinput/multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt" |>
#    read_tsv()
#
#fastqc_df <- 
#    left_join(meta_long, fastqc, join_by(fastq == Sample)) |>
#    pivot_longer(`Unique Reads`:`Duplicate Reads`, names_to = "read_type", values_to = "n") |>
#    mutate(read_type = sub("\\sReads$", "", read_type),
#           fastq = sub("^.+([12])$", "\\1", fastq),
#           id = paste(sample_id, stim, time, paste0("L", lane), paste0("fq", fastq), sep = "_"))
#
#keep_samples <- fastqc_df |>
#    filter(fastq == 1, sample_id != "BLANK") |>
#    group_by(sample_id, stim, time) |>
#    mutate(n = sum(n)) |>
#    ungroup() |>
#    filter(n >= 2e6) |>
#    distinct(plate, well, sample_id, stim, time) |>
#    unite("id", c(plate, well), sep = "_", remove = FALSE)
#
#keep_samples |>
#    mutate(sample_id = sub("\\.rep\\.\\d", "", sample_id)) |>
#    unite("sample_name", c("sample_id", "stim", "time"), sep = "_") |>
#    select(sample_id = id, sample_name) |>
#    write_tsv("./data/sample_decode.tsv")


# Salmon with replicates pooled at the fastq for alignment
meta_pooled <- 
    read_tsv("./data/metadata_pooledreps.tsv", col_names = c("sample_id", "fastq"))

salmon_pooled <- 
    sprintf("./results/salmon_pooledreps/%s/quant.sf", meta_pooled$sample_id) |>
    setNames(meta_pooled$sample_id) |>
    map_dfr(~read_tsv(.) |> 
	    left_join(gene_tx, join_by(Name == tx_id)),
	    .id = "sample_id") |>
    select(sample_id, gene_id, gene_name, tx_id = Name, read_count = NumReads, tpm = TPM)
	    
salmon_pooled_gene <- 
    salmon_pooled |>
    group_by(sample_id, gene_id, gene_name) |>
    summarise_at(vars(read_count, tpm), sum) |> 
    ungroup()

write_tsv(salmon_pooled, "./data/expression_pooled_reps_transcriptlevel.tsv")
write_tsv(salmon_pooled_gene, "./data/expression_pooled_reps.tsv")
