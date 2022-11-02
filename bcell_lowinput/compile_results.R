library(tidyverse)

# QC
meta_long <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_lowinput/metadata/description_hsapiens_longformat.tsv" %>%
    read_tsv() %>%
    pivot_longer(fq1:fq2, names_to = "dummy", values_to = "fastq") %>%
    mutate(fastq = basename(fastq),
           fastq = sub("\\.fastq\\.gz", "", fastq)) %>%
    select(-barcode_seq, -dummy)

fastqc_1 <- 
    "/home/ch229163/labshr/Lab_datasets/B_cells_lowinput/SK-56QR/get.broadinstitute.org/pkgs/SN0263576/multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt" %>%
    read_tsv()

fastqc_2 <- 
    "/home/ch229163/labshr/Lab_datasets/B_cells_lowinput/SK-56QS/get.broadinstitute.org/pkgs/SN0263542/multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt" %>%
    read_tsv()

fastqc_df <- bind_rows(fastqc_1, fastqc_2) %>%
    left_join(meta_long, ., c("fastq" = "Sample")) %>%
    pivot_longer(`Unique Reads`:`Duplicate Reads`, names_to = "read_type", values_to = "n") %>%
    mutate(read_type = sub("\\sReads$", "", read_type),
           fastq = sub("^.+([12])$", "\\1", fastq),
           id = paste(sample_id, stim, time, paste0("L", lane), paste0("fq", fastq), sep = "_"))

# Expression data
keep_samples <- fastqc_df %>%
    filter(fastq == 1, sample_id != "BLANK") %>%
    group_by(sample_id, stim, time) %>%
    mutate(n = sum(n)) %>%
    ungroup() %>%
    filter(n >= 2e6) %>%
    distinct(plate, well, sample_id, stim, time) 
    
meta_df <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_lowinput/metadata/description_hsapiens.tsv" %>%
    read_tsv() %>%
    select(plate:time) %>%
    inner_join(keep_samples) %>%
    unite("id", c(plate, well), sep = "_", remove = FALSE)

gene_tx <- "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.txToGene.tsv" %>%
    read_tsv()

quant_df <- file.path("./results/salmon", meta_df$id, "quant.sf") %>%
    setNames(meta_df$id) %>%
    map_df(read_tsv, .id = "id")

quant_gene <- quant_df %>%
    left_join(gene_tx, by = c("Name" = "transcript_id")) %>%
    group_by(id, gene_id, gene_name) %>%
    summarise(tpm = sum(TPM),
              count = sum(NumReads)) %>%
    ungroup()

write_rds(quant_gene, "./data/expression.rds")

meta_df %>%
    mutate(sample_id = sub("\\.rep\\.\\d", "", sample_id)) %>%
    unite("sample_name", c("sample_id", "stim", "time"), sep = "_") %>%
    select(sample_id = id, sample_name) %>%
    write_tsv("./data/sample_decode.tsv")
