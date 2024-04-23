library(tidyverse)
library(tximport)
library(glue)
library(patchwork)
library(ggrepel)
library(tidytext)

# Sample meta data
meta_data <- read_tsv("../0-qc_rnaseq/metadata.tsv")

stim_colors <- 
    c("Day 0" = "#898e9f",
      "BCR" = "#003967",
      "TLR7" = "#637b31",
      "DN2" = "#a82203")

# Transcript to Gene map
tx_to_gene <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf.gz") |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X3 == "transcript") |>
    mutate(tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	   gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(tx_id, gene_id, gene_name) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

# Expression data
salmon_files <-
    meta_data |>
    unite("prefix", c(donor_id, replic_id, stim), sep = "_") |>
    mutate(f = glue("./results/{prefix}/quant.sf")) |>
    select(prefix, f) |>
    deframe()

# Proportion of uniquely mapped reads
read_star_log <- function(f) {

    tmp <- 
	read_lines(f, skip = 7) |>
	trimws() |>
	{function(x) split(x, cumsum(grepl(":$", x)))}() |>
	{function(x) setNames(x, map_chr(x, 1))}() |>
	map(~.x[-1]) |>
	map_dfr(~tibble(x = .) |> 
		separate(x, c("info", "value"), sep = "\\|\t"), 
		.id = "reads") |>
        mutate(reads = str_remove(reads, " READS:"),
	       info = trimws(info),
	       value = parse_number(value))

    tmp |>
	filter(grepl("reads number|^Number of (chimeric )?reads", info)) |>
	mutate(value = as.integer(value)) |>
	select(reads, info, value) |>
	group_by(reads) |>
	summarise(number_of_reads = sum(value)) |>
	ungroup() |> 
	mutate(reads = tolower(reads)) |>
	arrange(reads)
}

temp_work <- system("echo $TEMP_WORK", intern = TRUE)

star_log_df <- 
    meta_data |>
    unite("sample_id", c(donor_id, replic_id), sep = "_") |>
    mutate(logfile = file.path(temp_work, glue("/bam/highinput_merged/{sample_id}_{stim}_Log.final.out")),
	   data = map(logfile, read_star_log)) |>
    mutate(stim = recode(stim, "unstday0" = "Day 0"),
	   stim = factor(stim, levels = names(stim_colors))) |>
    select(sample_id, stim, data) |>
    arrange(sample_id, stim) |>
    unnest(cols = data) |>
    group_by(reads) |>
    filter(any(number_of_reads) > 0) |>
    mutate(m = mean(number_of_reads)) |>
    ungroup() |>
    arrange(sample_id, stim, desc(m)) |>
    mutate(reads = fct_inorder(reads)) |>
    select(-m)

perc_uniq_plot <- 
    ggplot(star_log_df,
	   aes(x = number_of_reads, y = sample_id)) +
    geom_col(aes(fill = stim, alpha = fct_rev(reads)), 
	     position = "fill", color = "black", linewidth = .2) +
    scale_x_continuous(breaks = c(0, .5, 1)) +
    scale_fill_manual(values = stim_colors) +
    scale_alpha_manual(values = rev(c(1, .6, .2))) +
    facet_wrap(~stim, ncol = 2, scales = "free_y") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = "Percentage of reads", 
	 y = NULL,
	 alpha = "Reads:") +
    guides(fill = "none")

ggsave("./plots/perc_uniq.png", perc_uniq_plot)


# TPMs from Salmon
salmon_meta <- meta_data |>
    select(1:3) |>
    unite("sample_id", c(donor_id, replic_id, stim), sep = "_", remove = FALSE) |>
    unite("sample_id2", c(donor_id, replic_id), sep = "_", remove = FALSE) |>
    select(sample_id, sample_id2, stim)

salmon_df <-
    salmon_files |>
    map_dfr(~. |>
	    read_tsv() |>
	    left_join(tx_to_gene, join_by(Name == "tx_id")) |>
	    group_by(gene_id, gene_name) |>
	    summarise(tpm = sum(TPM)) |>
	    ungroup(),
	    .id = "sample_id") |>
    left_join(salmon_meta, join_by(sample_id)) |>
    select(sample_id = sample_id2, stim, gene_id, gene_name, tpm)

const_expr_genes <- 
    salmon_df |>
    group_by(gene_id) |>
    filter(mean(tpm >= 1) >= .9) |>
    ungroup() |>
    distinct(gene_id, gene_name)

qc_plot_data <- 
    salmon_df |>
    inner_join(const_expr_genes, join_by(gene_id, gene_name)) |>
    group_by(sample_id, stim) |>
    summarise(m = mean(tpm >= 1)) |>
    ungroup() |>
    arrange(m) |>
    mutate(stim = recode(stim, "unstday0" = "Day 0"),
	   stim = factor(stim, levels = names(stim_colors)))
   
qc_plot <- 
    ggplot(qc_plot_data, aes(x = m, y = sample_id)) +
    geom_col(aes(fill = stim), 
	     color = "black", linewidth = .2) +
    geom_vline(xintercept = c(.96, .98, 1), linetype = 2, linewidth = .5) +
    scale_x_continuous(limits = c(0, 1),
		       breaks = c(.96, .98, 1),
		       labels = c("0.96", "0.98", "1")) +
    coord_cartesian(xlim = c(0.96, 1)) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 2, scales = "free_y") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  strip.text = element_text(size = 11),
	  axis.text.x = element_text(size = 10),
	  axis.title.x = element_text(size = 10),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = "Proportion of genes expressed", 
	 y = NULL,
	 title = "Proportion of usually expressed genes*\nthat are expressed by each sample",
	 subtitle = "* Genes with at least 1 TPM in 90% of samples") +
    guides(fill = "none")
    
ggsave("./plots/const_genes.png", qc_plot, width = 5, height = 7.5)

# PCA with DESeq2
library(DESeq2)

# tximport method
txi <- 
    tximport(salmon_files, 
	     type = "salmon", 
	     tx2gene = select(tx_to_gene, tx_id, gene_id))

sample_table <- 
    tibble(sample_id = colnames(txi$counts)) |>
    extract(sample_id, c("group"), "([^_]+)$", remove = FALSE) |>
    column_to_rownames("sample_id")

dds <- 
    DESeqDataSetFromTximport(txi, sample_table, ~1) |>
    estimateSizeFactors() |>
    {function(x) x[rowSums(counts(x, normalized = TRUE) >= 5 ) >= 18, ]}() |>
    vst()

top_variable_genes <- 
    rowVars(assay(dds), useNames = TRUE) |> 
    enframe("gene_id", "var") |>
    arrange(desc(var)) |>
    dplyr::slice(1:5000)

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

pca_df <- 
    pc_scores |>
    select(sample_id, PC1:PC4) |>
    extract(sample_id, c("sample_id", "stim"), "(\\d+_[123])_(.+)") |>
    mutate(stim = recode(stim, "unstday0" = "Day 0"),
	   stim = factor(stim, levels = names(stim_colors)))

labels_samples <- 
    tribble(~sample_id, ~stim,
	    "10098629_1", "DN2",
	    "10095189_1", "DN2",
	    "10098629_1", "TLR7",
	    "10095189_1", "TLR7",
	    "10050385_1", "DN2",
	    "10050385_1", "BCR",
	    "10050385_1", "TLR7",
	    "10068703_1", "Day 0",
	    "10028815_1", "Day 0",
	    "10098629_1", "Day 0")


pca_vst_1 <- 
    ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(color = stim), size = 2) +
    geom_label_repel(data = inner_join(pca_df, labels_samples),
		    aes(label = sample_id, color = stim), 
		    size = 2.5,
		    min.segment.length = .2,
		    segment.size = .2,
		    max.overlaps = Inf,
		    alpha = .5,
		    label.padding = 0.1,
		    show.legend = FALSE) +
    scale_color_manual("Stim:", values = stim_colors) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white"),
	  legend.position.inside = c(.8, .8),
	  legend.box.background = element_rect(color = "black", linewidth = .1)) +
    labs(x = pc_varexp$lab[1], y = pc_varexp$lab[2]) +
    guides(color = guide_legend(position = "inside"))

ggsave("./plots/pca_vst.png", pca_vst_1, height = 4, width = 4)


