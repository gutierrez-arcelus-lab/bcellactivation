library(tidyverse)
library(tximport)
library(DESeq2)

slice <- dplyr::slice

# Sample meta data
meta_data <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/metadata.tsv" |>
    read_tsv()

# Transcript to Gene map
tx_to_gene <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf") |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X3 == "transcript") |>
    mutate(tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	   gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(tx_id, gene_id, gene_name)

# Expression data
salmon_files <- 
    sprintf("./results/salmon/%s_%s/quant.sf", meta_data$sample_id, meta_data$stim) |>
    setNames(paste(meta_data$sample_id, meta_data$stim, sep = "_"))

txi <- tximport(salmon_files, 
		type = "salmon", 
		tx2gene = select(tx_to_gene, tx_id, gene_id))

sample_table <- 
    tibble(sample_id = colnames(txi$counts)) |>
    extract(sample_id, c("group"), "[^_]+_(.+)", remove = FALSE) |>
    column_to_rownames("sample_id")

dds <- 
    DESeqDataSetFromTximport(txi, sample_table, ~1) |>
    estimateSizeFactors() |>
    {function(x) x[rowSums(counts(x, normalized = TRUE) >= 10 ) >= 3, ]}() |>
    vst()

top_variable_genes <- 
    rowVars(assay(dds)) |> 
    setNames(rownames(assay(dds))) |>
    enframe("gene_id", "var") |>
    arrange(desc(var)) |>
    slice(1:5000)

count_matrix <- t(assay(dds)[top_variable_genes$gene_id, ])

pca <- prcomp(count_matrix, center = TRUE, scale. = TRUE)

pc_scores <- as_tibble(pca$x, rownames = "sample_id")

pc_eigenvals <- pca$sdev^2

pc_varexp <- 
    tibble(PC = paste0("PC", 1:length(pc_eigenvals)),
	   variance = pc_eigenvals) |>
    mutate(pct = round(variance/sum(variance) * 100, 1),
           pct = paste0("(", pct, "%)")) |>
    select(PC, pct) |>
    unite("lab", c(PC, pct), sep = " ", remove = FALSE)

stim_colors <- 
    c("Day 0" = "grey30",
      "BCR" = "#4DBBD5FF",
      "TLR7" = "#91D1C2FF",
      "DN2" = "#E64B35FF")

match_mbv <- 
    read_tsv("./3-mbv/match_summary.tsv", col_types = c(.default = "c")) |>
    mutate(m = bam_donor_id == vcf_donor_id) |>
    select(sample_id, stim = bam_stim, m)

pca_df <- 
    pc_scores |>
    select(sample_id, PC1:PC2) |>
    separate(sample_id, c("sample_id", "stim"), sep = "_") |>
    left_join(match_mbv) |>
    mutate(stim = recode(stim, "unstday0" = "Day 0"),
	   stim = fct_inorder(stim))

pca_vst <- 
    ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(data = filter(pca_df, m),
	       aes(fill = stim), size = 4, shape = 21, alpha = .25) +
    geom_point(data = filter(pca_df, !m),
	       aes(fill = stim), size = 4, shape = 21, alpha = 1) +
    scale_fill_manual("Condition:", values = stim_colors) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.key.height = unit(.25, "cm"),
	  legend.text = element_text(size = 8)) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = pc_varexp$lab[1], y = pc_varexp$lab[2])

ggsave("./plots/pca_vst.png", pca_vst, height = 5.25)

