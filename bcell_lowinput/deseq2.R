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
    {function(x) x[rowSums(counts(x, normalized = TRUE) >= 1 ) >= 3, ]}() |>
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

