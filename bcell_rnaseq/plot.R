library(tidyverse)
library(tidytext)
library(broom)
library(cowplot)

stim_colors <- c("unstday0" = "grey",
		 "BCR" = "#4DBBD5FF",
		 "TLR7" = "#91D1C2FF",
		 "DN2" = "#E64B35FF")
# QC

qc_df <- 
    "/home/ch229163/labshr/Lab_datasets/B_cells_rnaseq/QC/multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt" |>
    read_tsv() |>
    extract(Sample, 
	    c("subject_id", "stim", "r"), 
	    "\\d+_(\\d+[_123]*)_(.+)_MG\\d+_.+_(R[12])") |>
    separate(subject_id, c("subject_id", "sample_id"), sep = "_") |>
    mutate(sample_id = replace_na(sample_id, "1"),
	   sample_id = paste(subject_id, sample_id, sep = ".")) |>
    pivot_longer(`Unique Reads`:`Duplicate Reads`, names_to = "read_type", values_to = "n") |>
    mutate(read_type = sub("\\sReads$", "", read_type))

plot_df <- qc_df %>%
    filter(r == "R1") %>%
    mutate(stim = recode(stim, "TR7" = "TLR7"),
	   stim = factor(stim, levels = names(stim_colors)))

total_reads <- ggplot(plot_df, 
       aes(x = reorder_within(sample_id, by = n, within = stim, fun = min), 
	   y = n)) +
    geom_col(aes(fill = stim, alpha = read_type),
	     color = "black", size = .1, width = 1) +
    scale_fill_manual(values = stim_colors) +
    scale_alpha_manual(values = c("Unique" = 1, "Duplicate" = .5)) +
    scale_x_reordered() +
    scale_y_continuous(labels = scales::comma, breaks = scales::pretty_breaks(6)) +
    facet_wrap(~stim, scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey80", size = .5),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
	  axis.text.y = element_text(size = 14),
	  axis.title = element_text(size = 16),
	  strip.text = element_text(size = 16),
          plot.background = element_rect(fill = "white", color = "white"),
          legend.position = "none") +
    labs(x = NULL, y = "Total reads", fill = "Stim")

ggsave("./plots/total_reads.png", total_reads, width = 10, height = 5)



# PCA 
quant <- read_tsv("./quantifications_genes.tsv")

variable_genes <- quant |>
    group_by(gene_id, gene_name) |>
    summarise(v = var(tpm)) |>
    ungroup() |>
    top_n(2000, v) |>
    arrange(desc(v))

variable_matrix <- quant |>
    filter(gene_id %in% variable_genes$gene_id) |>
    select(-gene_name, -counts) |>
    pivot_wider(names_from = gene_id, values_from = tpm) |>
    unite("id", c(sample_id, stim), sep = "_") |>
    column_to_rownames("id")

pca <- prcomp(variable_matrix, center = TRUE, scale. = TRUE, rank. = 10)

var_exp <- pca |>
    tidy("pcs") |>
    filter(PC %in% 1:10) |>
    ggplot(aes(x = factor(PC), y = percent)) +
	geom_col(alpha = 0.7) +
	scale_y_continuous(labels = scales::label_percent()) +
    theme_minimal() +
    theme(text = element_text(size = 14),
	  panel.grid.major.x = element_blank()) +
    labs(x = "PC", y = "Variance explained")

pc_scores <- as_tibble(pca$x, rownames = "id")

pca_df <- pc_scores %>%
    separate(id, c("sample_id", "stim"), sep = "_") %>%
    mutate(stim = factor(stim, levels = names(stim_colors)))

pca_plot <- ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(fill = stim), size = 4, shape = 21) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(text = element_text(size = 14),
	  panel.grid = element_blank(),
          legend.key.height = unit(.25, "cm")) +
    guides(fill = guide_legend(ncol = 1))

pca_out <- plot_grid(var_exp, pca_plot, nrow = 1, rel_widths = c(.5, 1)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/pca.png", pca_out, width = 8, height = 5)


# Edger
cpm_df <- read_tsv("./edger_cpm.tsv") |>
    mutate(stim_test = factor(stim_test, levels = c("BCR", "TLR7", "DN2")),
	   stim = factor(stim, levels = c("unstday0", "BCR", "TLR7", "DN2")))

edger_res <- read_tsv("./edger_results.tsv") |>
    mutate(stim_test = factor(stim_test, levels = c("BCR", "TLR7", "DN2")))

sle_genes <- read_tsv("../bcell_scrna/reported_genes.tsv")

plot_fc <- edger_res |> 
    filter(gene_name %in% sle_genes$gene, FDR < 0.05) |>
    group_by(stim_test) |>
    top_n(25, abs(logFC)) |>
    ungroup() |>
    ggplot(aes(x = logFC,y = reorder_within(gene_name, within = stim_test, by = logFC))) +
	geom_col(aes(fill = logFC)) +
	scale_y_reordered() +
	scale_fill_gradient2() +
	facet_wrap(~stim_test, nrow = 1, scales = "free") +
	theme_bw() +
	theme(legend.position = "none",
	      panel.grid.minor.x = element_blank()) +
	labs(y = NULL)

ggsave("./plots/fold_change.png", plot_fc, width = 10, height = 5)


stat1 <- cpm_df |> 
    filter(gene_name == "STAT1") |>
    ggplot(aes(x = stim, y = cpm)) +
	geom_line(aes(group = sample_id), linewidth = .2) +
	geom_point(aes(fill = stim), size = 4, shape = 21, stroke = .2) +
	scale_fill_manual(values = stim_colors) +
	facet_wrap(~stim_test, scales = "free_x", nrow = 1) +
	theme_bw() +
	theme(legend.position = "none",
	      panel.grid.minor.y = element_blank()) +
	labs(x = NULL, y = "Counts per Million")

ggsave("./plots/stat1.png", stat1, height = 3)


irf8 <- cpm_df |> 
    filter(gene_name == "IRF8") |>
    ggplot(aes(x = stim, y = cpm)) +
	geom_line(aes(group = sample_id), linewidth = .2) +
	geom_point(aes(fill = stim), size = 4, shape = 21, stroke = .2) +
	scale_fill_manual(values = stim_colors) +
	facet_wrap(~stim_test, scales = "free_x", nrow = 1) +
	theme_bw() +
	theme(legend.position = "none",
	      panel.grid.minor.y = element_blank()) +
	labs(x = NULL, y = "Counts per Million")

ggsave("./plots/irf8.png", irf8, height = 3)

# GO

go_res <- read_tsv("./top_diff_go.tsv") |>
    mutate(stim = factor(stim, levels = names(stim_colors[-1]))) |>
    arrange(desc(-log10(P.Up))) |>
    mutate(Term = fct_inorder(Term))

goplot <- ggplot(go_res) +
    geom_hline(yintercept = -log10(0.05), linewidth = .5) +
    geom_hline(yintercept = -log10(0.01), linewidth = .5, linetype = "dashed") +
    geom_col(aes(x = stim, y = -log10(P.Up), fill = stim), position = "dodge") +
    scale_fill_manual(values = stim_colors[-1]) +
    facet_wrap(~Term, nrow = 4,
	       labeller = labeller(Term = label_wrap_gen(width = 25))) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = NULL)

ggsave("./plots/go.png", goplot, width = 10, height = 6)
