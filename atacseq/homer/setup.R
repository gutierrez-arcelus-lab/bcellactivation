library(tidyverse)
library(glue)
     
# select negative fold-change because analysis was done using unstim as effect stim
process_bed <- function(s) {

    glue("../results_deseq2/{s}_24vsunst_24.tsv") |>
    read_tsv() |>
    filter(!is.na(padj), padj <= 0.01, log2FoldChange >= 1) |>
    select(Chr, Start, End, interval, log2FoldChange, Strand) |>
    mutate(log2FoldChange = round(log2FoldChange, 2)) |>
    write_tsv(glue("./data/{s}.bed"), col_names = FALSE)
}

stims <- read_lines("./stims.txt")

walk(stims, process_bed)

# Use consensus peaks for background, instead of random genomics sequences
consensus <- 
    "../results/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.bed" |>
    read_tsv(col_names = FALSE) |>
    arrange(X1, X2)

write_tsv(consensus, "./data/consensus.bed", col_names = FALSE)

