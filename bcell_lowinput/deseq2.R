library(tidyverse)
library(DESeq2)
library(BiocParallel)
library(furrr)

cpus <- availableCores() |> as.integer()

gene_names <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf.gz") |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X3 == "gene") |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(gene_id, gene_name)

dat <- read_rds("./wgcna/data/gene_expression.rds")

sample_table <- 
    tibble(sample_id = colnames(dat$counts)) |>
    separate(sample_id, c("donor_id", "stim", "timept"), sep = "_", remove = FALSE) |>
    unite("group", c(stim, timept), sep = "_") |>
    column_to_rownames("sample_id")

dds <- 
    DESeqDataSetFromTximport(dat, sample_table, ~group + donor_id) |>
    estimateSizeFactors() |>
    DESeq(parallel = TRUE, BPPARAM = MulticoreParam(cpus))



comp0 <- results(dds, contrast = c("group", "DN2_72hrs", "Unstim_0hrs"), alpha = 0.05)

comp0 |> 
    as_tibble(rownames = "gene_id") |> 
    left_join(gene_names) |>
    filter(padj <= 0.05) |>
    arrange(padj)






comp1 <- results(dds, contrast = c("group", "BCR_72hrs", "DN2_72hrs"), alpha = 0.05)

comp1 |> 
    as_tibble(rownames = "gene_id") |> 
    left_join(gene_names) |>
    filter(padj <= 0.05, log2FoldChange < 0) |>
    arrange(padj)


comp2 <- results(dds, contrast = c("group", "BCR_4hrs", "DN2_4hrs"), alpha = 0.05)

comp2 |> 
    as_tibble(rownames = "gene_id") |> 
    left_join(gene_names) |>
    filter(padj <= 0.05) |>
    arrange(padj)




# Save counts
sample_info <- 
    tibble(sample_id = colnames(counts(dds))) |>
    separate(sample_id, 
	     c("donor", "stim", "timep"), 
	     sep = "_", remove = FALSE) |>
    mutate(timep = parse_number(timep),
	   timep = factor(timep, levels = sort(unique(timep))),
	   stim = recode(stim, "IL4" = "IL-4c", "CD40L" = "CD40c", 
			 "TLR9" = "TLR9c", "TLR7" = "TLR7c",
			 "BCR" = "BCRc", "BCR-TLR7" = "BCR/TLR7c", "DN2" = "DN2c"),
	   stim = factor(stim, levels = c("Unstim", "IL-4c", "CD40c", "TLR7c", 
					  "TLR9c", "BCRc", "BCR/TLR7c", "DN2c"))) |>
    arrange(stim, timep, donor)

gene_ids <-
    gene_names |> 
    add_count(gene_name) |>
    mutate(gene_label = case_when(n == 1 ~ gene_name,
				  n > 1 ~ paste(str_remove(gene_id, "\\.\\d+$"), gene_name, sep = "_"))) |>
    select(gene_id, gene_label)


counts_df <- 
    counts(dds, normalized = TRUE) |>
    as.data.frame() |>
    rownames_to_column("gene_id") |>
    as_tibble() |>
    left_join(gene_ids, join_by(gene_id)) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    select(gene_id, gene_label, everything()) |>
    pivot_longer(-(gene_id:gene_label), names_to = "sample_id", values_to = "norm_counts") |>
    left_join(sample_info, join_by(sample_id)) |>
    select(donor, stim, timep, gene_label, norm_counts) |>
    arrange(stim, timep, donor, gene_label)



write_rds(counts_df, "./data/deseq_normalized_counts.rds")
#

