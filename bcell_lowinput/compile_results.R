library(tidyverse)

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

fastqc <-
    "/temp_work/ch229163/fastq/lowinput/multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt" |>
    read_tsv()

fastqc_df <- 
    left_join(meta_long, fastqc, join_by(fastq == Sample)) |>
    pivot_longer(`Unique Reads`:`Duplicate Reads`, names_to = "read_type", values_to = "n") |>
    mutate(read_type = sub("\\sReads$", "", read_type),
           fastq = sub("^.+([12])$", "\\1", fastq),
           id = paste(sample_id, stim, time, paste0("L", lane), paste0("fq", fastq), sep = "_"))

# Expression data
keep_samples <- fastqc_df |>
    filter(fastq == 1, sample_id != "BLANK") |>
    group_by(sample_id, stim, time) |>
    mutate(n = sum(n)) |>
    ungroup() |>
    filter(n >= 2e6) |>
    distinct(plate, well, sample_id, stim, time) |>
    unite("id", c(plate, well), sep = "_", remove = FALSE)
    
gene_tx <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf.gz" |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X3 == "transcript") |>
    mutate(tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	   gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(tx_id, gene_id, gene_name)

quant_df <- file.path("./results/salmon", keep_samples$id, "quant.sf") |>
    setNames(keep_samples$id) |>
    map_df(read_tsv, .id = "id")

quant_tx <- quant_df |>
    left_join(gene_tx, by = c("Name" = "tx_id")) |>
    select(id, gene_id, gene_name, tx_id = Name, count = NumReads, tpm = TPM) |>
    group_by(tx_id) |>
    filter(!all(tpm == 0)) |>
    ungroup()

write_rds(quant_tx, "./data/expression_transcripts.rds")

quant_gene <- quant_tx |>
    group_by(id, gene_id, gene_name) |>
    summarise_at(vars(count, tpm), sum) |> 
    ungroup()

write_rds(quant_gene, "./data/expression.rds")

keep_samples |>
    mutate(sample_id = sub("\\.rep\\.\\d", "", sample_id)) |>
    unite("sample_name", c("sample_id", "stim", "time"), sep = "_") |>
    select(sample_id = id, sample_name) |>
    write_tsv("./data/sample_decode.tsv")

# Salmon with replicates pooled at the fastq for alignment
meta_pooled <- 
    read_tsv("./data/metadata_pooledreps.tsv", col_names = c("sample_id", "fastq"))

salmon_pooled <- 
    sprintf("./results/salmon_pooledreps/%s/quant.sf", meta_pooled$sample_id) |>
    setNames(meta_pooled$sample_id) |>
    map_dfr(~read_tsv(.) |> 
	    left_join(gene_tx, join_by(Name == tx_id)) |>
	    group_by(gene_id, gene_name) |>
	    summarise(count = sum(NumReads),
		      tpm = sum(TPM)) |>
	    ungroup(),
	    .id = "sample_id") |>
    separate(sample_id, c("sample_id", "stim", "hours"), sep = "_")

write_tsv(salmon_pooled, "./data/expression_pooled_reps.tsv")

