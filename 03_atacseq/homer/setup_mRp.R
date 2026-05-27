library(tidyverse)
library(glue)

dir.create("./data_mRp")

process_bed <- function(s) {

    glue("../results_deseq2/mRp/{s}_24vsunst_24.tsv") |>
    read_tsv() |>
    filter(!is.na(padj), padj <= 0.01, log2FoldChange >= 1) |>
    select(Chr, Start, End, interval, log2FoldChange, Strand) |>
    mutate(log2FoldChange = round(log2FoldChange, 2)) |>
    write_tsv(glue("./data_mRp/{s}.bed"), col_names = FALSE)
}

stims <- read_lines("./stims.txt")

walk(stims, process_bed)

# Use consensus peaks for background, instead of random genomics sequences
consensus <- 
    "../results/bwa/merged_replicate/macs2/narrow_peak/consensus/consensus_peaks.mRp.clN.bed" |>
    read_tsv(col_names = FALSE) |>
    arrange(X1, X2)

write_tsv(consensus, "./data_mRp/consensus.bed", col_names = FALSE)

