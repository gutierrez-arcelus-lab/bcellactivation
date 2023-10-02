library(tidyverse)

sample_specs <- read_csv("./samplesheet.csv") |>
    select(stim = sample, donor_id = fastq_1, replicate) |>
    mutate(donor_id = basename(donor_id)) |>
    extract(donor_id, "donor_id", "\\d+_([^_]+)_.+") |>
    distinct() |>
    filter(donor_id != "3donors", stim != "IL4_24") |>
    mutate(sample_id = paste0(stim, "_REP", replicate))

count_table <- 
    "./results/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.featureCounts.txt" |>
    read_tsv(skip = 1)

colnames(count_table) <- gsub("\\.mLb\\.\\clN\\.sorted\\.bam", "", colnames(count_table))

peak_annot <- 
    "./results/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.annotatePeaks.txt" |> 
    read_tsv() |>
    select(peak_id = 1, everything())



atac <- count_table |>
    select(peak_id = Geneid, !!!sample_specs$sample_id) |>
    pivot_longer(-peak_id, names_to = "sample_id") |>
    group_by(sample_id) |>
    mutate(rppm = (value/sum(value)) * 1e6) |>
    group_by(peak_id) |> 
    filter(!all(rppm == 0)) |>
    ungroup() |>
    select(peak_id, sample_id, rppm)

# RNA-seq
# use only rep "1" for samples with multiple technical replicates
# An alternative would be to merge replicates at the aligment level such as the ATAC-seq pipeline
rna <- read_tsv("../bcell_rnaseq/results/salmon_genes.tsv", col_types = "ccccdd") |>
    filter(grepl("\\.1$", sample_id)) |>
    mutate(sample_id = sub("\\.1$", "", sample_id),
	   stim = recode(stim, "unstday0" = "unst_0", "BCR" = "BCR_24", 
			 "TLR7" = "TLR7_24", "DN2" = "DN2_24")) |>
    select(donor_id = sample_id, stim, gene_id, gene_name, tpm) |>
    inner_join(sample_specs, join_by(donor_id, stim)) |> 
    group_by(gene_id) |> 
    filter(!all(tpm == 0)) |>
    ungroup()


cor_df <- 
    rna |>
    filter(gene_name == "IKZF1") |>
    left_join(atac, join_by(sample_id), multiple = "all") |>
    select(donor_id, sample_id, stim, gene_id, gene_name, peak_id, tpm, rppm) |>
    group_by(gene_id, gene_name, peak_id) |>
    reframe(r = as_tibble(cor.test(tpm, rppm)[c("estimate", "p.value")])) |>
    unnest(cols = r)

cor_df |> filter(p.value < 0.05/nrow(cor_df))
cor_df |> filter(p.adjust(p.value) < 0.05)
cor_df |> arrange(desc())

peak_annot |> 
    filter(peak_id == "Interval_93698") |> 
    print(width = Inf)

