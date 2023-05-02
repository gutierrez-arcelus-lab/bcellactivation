library(tidyverse)

sample_specs <- read_csv("./samplesheet.csv") |>
    select(stim = sample, donor_id = fastq_1, replicate) |>
    mutate(donor_id = basename(donor_id)) |>
    extract(donor_id, "donor_id", "\\d+_([^_]+)_.+") |>
    distinct() |>
    filter(donor_id != "3donors") |>
    mutate(sample_id = paste0(stim, "_REP", replicate))

count_table <- 
    "./results/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.featureCounts.txt" |>
    read_tsv(skip = 1)

colnames(count_table) <- gsub("\\.mLb\\.\\clN\\.sorted\\.bam", "", colnames(count_table))

atac <- count_table |>
    select(peak_id = Geneid, !!!sample_specs$sample_id) |>
    pivot_longer(-peak_id, names_to = "sample_id") |>
    group_by(sample_id) |>
    mutate(rppm = (value/sum(value)) * 1e6) |>
    ungroup() |>
    select(peak_id, sample_id, rppm)

rna <- read_tsv("../bcell_rnaseq/results/salmon_genes.tsv", col_types = "ccccdd")

# combine replicates
# recompute TPMs?

rna |> 
    distinct(sample_id, stim) |>
    separate(sample_id, c("donor_id", "rep"), sep = "\\.") |>
    filter(donor_id %in% sample_specs$donor_id)
