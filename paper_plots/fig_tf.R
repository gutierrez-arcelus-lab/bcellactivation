library(tidyverse)
library(patchwork)

gtf_file <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf.gz"

gtf_import <-
    rtracklayer::import(gtf_file)

tx2gene <-
    gtf_import |>
    as_tibble() |>
    filter(type == "transcript") |>
    select(tx_id = transcript_id, gene_id, gene_name) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

meta_data <-
    "../bcell_lowinput/data/metadata_pooledreps.tsv" |>
    read_tsv(col_names = c("sample_id", "fq1", "fq2")) |>
    separate(sample_id, c("donor", "stim", "timep"), sep = "_", remove = FALSE) |>
    select(sample_id, donor, stim, timep)

salmon_files <-
    file.path("../bcell_lowinput/results/salmon_pooledreps", meta_data$sample_id, "quant.sf") |>
    setNames(meta_data$sample_id)

# Salmon
salmon_df <- 
    map_dfr(salmon_files, 
	    ~read_tsv(.) |> 
	    left_join(tx2gene, join_by(Name == tx_id)) |>
	    group_by(gene_id, gene_name) |>
	    summarize(tpm = sum(TPM),
		      gene_count = sum(NumReads)
	    ) |>
	    ungroup(),
	    .id = "sample_id")

library_sizes <- 
    salmon_df |> 
    group_by(sample_id) |>
    summarize(total_reads = sum(gene_count)) |>
    ungroup() |>
    arrange(total_reads)

samples_pass <-
    filter(library_sizes, total_reads >= 2e6) |>
    pull(sample_id)

# TFs
tf_df <- 
    "https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.txt" |>
    read_tsv() |>
    janitor::clean_names() |> 
    filter(is_tf == "Yes")

diff_expr <- 
    "../bcell_lowinput/results/edger/diff_expr_all_times.tsv" |>
    read_tsv()

deg_tfs <- 
    diff_expr |>
    filter(group1 %in% c("Unstim.4", "Unstim.24"), 
	   group2 %in% c("IL4.4", "IL4.24", "TLR7.4", "TLR7.24", "BCR.4", "BCR.24", "DN2.4", "DN2.24")) |>
    separate(group1, c("stim1", "timep1"), sep = "\\.") |>
    separate(group2, c("stim2", "timep2"), sep = "\\.") |>
    filter(timep1 == timep2, logFC >= 0.5, FDR <= 0.05) |>
    distinct(gene_id) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    filter(gene_id %in% tf_df$ensembl_id) |>
    pull(gene_id)

#expressed_tfs <- 
#    salmon_df |>
#    filter(sample_id %in% samples_pass,
#	   gene_id %in% tf_df$ensembl_id) |>
#    left_join(meta_data, join_by(sample_id)) |>
#    group_by(stim, timep, gene_id, gene_name) |>
#    summarize(mean_tpm = mean(tpm)) |>
#    ungroup() |>
#    group_by(gene_id, gene_name) |>
#    filter(any(mean_tpm >= 10)) |>
#    ungroup() |>
#    distinct(gene_id, gene_name)

# DESeq2
txi <- tximport::tximport(salmon_files[samples_pass], type = "salmon", tx2gene = tx2gene)

sample_table <- 
    tibble(sample_id = colnames(txi$counts)) |>
    left_join(meta_data) |>
    column_to_rownames("sample_id")

dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ 1)
vsd <- DESeq2::vst(dds, blind = TRUE)
vst_mat <- SummarizedExperiment::assay(vsd)

vst_mat_tf <- vst_mat[rownames(vst_mat) %in% deg_tfs, ]

gene_dist <- dist(vst_mat_tf, method = "euclidean")
gene_hclust <- hclust(gene_dist, method = "complete")
gene_cluster_order <- rownames(vst_mat_tf)[gene_hclust$order]

sample_dist <- dist(t(vst_mat_tf), method = "euclidean")
sample_hclust <- hclust(sample_dist, method = "complete")
sample_cluster_order <- rownames(t(vst_mat_tf))[sample_hclust$order]

vst_mat_tf_scaled <- t(scale(t(vst_mat_tf)))

#vst_mat_tf_scaled2 <- 
#    t(vst_mat_tf) |>
#    as.data.frame() |>
#    mutate_all(GenABEL::rntransform) |>
#    t()

vst_tf_data <- 
    vst_mat_tf_scaled |>
    as.data.frame() |>
    rownames_to_column("gene_id") |>
    as_tibble() |>
    left_join(distinct(tx2gene, gene_id, gene_name), join_by(gene_id)) |>
    select(gene_id, gene_name, everything()) |>
    pivot_longer(-c(gene_id, gene_name), names_to = "sample_id") |>
    left_join(meta_data, join_by(sample_id)) |>
    mutate(timep = str_remove(timep, "hrs"),
	   timep = as.integer(timep),
	   stim = case_match(stim, 
			     "IL4" ~"IL-4c",
			     "CD40L" ~ "CD40c",
			     "TLR9" ~ "TLR9c",
			     "TLR7" ~ "TLR7c",
			     "BCR" ~ "BCRc",
			     "BCR-TLR7" ~ "BCR/TLR7c",
			     "DN2" ~ "DN2c",
			     .default = stim),
	   stim = factor(stim, levels = c("Unstim", "IL-4c", "CD40c", "TLR9c", "TLR7c", "BCRc", "BCR/TLR7c", "DN2c")),
	   gene_id = factor(gene_id, levels = gene_cluster_order),
	   sample_id = factor(sample_id, levels = sample_cluster_order)) |>
    arrange(gene_id, stim, timep) |>
    select(gene_id, gene_name, sample_id, stim, timep, value) |>
    mutate(gene_name = fct_inorder(gene_name))

heatplot <- 
    ggplot(vst_tf_data) +
    geom_tile(aes(x = sample_id, y = gene_name, fill = value)) +
    scale_fill_gradient2(
			 low = "#542788",  # Deep Purple
			 mid = "#F7F7F7", 
			 high = "#E08214", # Deep Orange
			 midpoint = 0,
			 limits = c(-3, 3),
			 oob = scales::squish
    ) +
    #scale_fill_gradient2(
    #    		 low = "#2166AC",    # blue
    #                     mid = "white",
    #                     high = "#B2182B",   # red
    #                     midpoint = 0
    #) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 4),
	  axis.text.x = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")
    ) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colorbar(barheight = 12, barwidth = .5))

# stim panel
stim_colors <- 
    read_tsv("./figure_colors.txt", col_names = c("stim", "timep", "color")) |>
    unite("cond", c(stim, timep), sep = " ") |>
    deframe()

stim_panel_df <- 
    vst_tf_data |>
    distinct(sample_id, stim, timep) |>
    arrange(sample_id) |>
    unite("cond", c(stim, timep), sep = " ")

stim_panel_plot <- 
    ggplot(stim_panel_df) +
    geom_tile(aes(x = sample_id, y = factor(1), fill = cond)) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  panel.grid = element_blank(),
	  plot.margin = margin(0, 0, 0, 0),
	  plot.background = element_rect(fill = "white", color = "white")
    ) +
    guides(fill = "none") +
    coord_cartesian(ylim = c(1, 1))


ggsave("./fig_tf.png", stim_panel_plot + heatplot + plot_layout(ncol = 1, heights = c(.01, 1)), width = 12, height = 18)
