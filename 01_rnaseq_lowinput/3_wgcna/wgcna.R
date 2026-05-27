# ==============================================================================
# Description:  Constructs weighted gene co-expression networks for specific 
#               stimulations. Filters for dynamically expressed genes (from edgeR), 
#               applies DESeq2 Variance Stabilizing Transformation (VST) to, and 
#               identifies co-expressed modules. 
#               Finally, performs GO enrichment on modules.
# Input:        1. gene_expression_txi.rds
#               2. results.tsv (edgeR time-course results)
#               3. Command line argument for specific stimulation (e.g., "DN2")
# Output:       Module assignments, eigengenes, kME values, trend plots, 
#               and Gene Ontology enrichment results.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment & Setup
# ------------------------------------------------------------------------------
library(tidyverse)
library(glue)
library(furrr)
library(tximport)
library(DESeq2)
library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(tidytext)
library(patchwork)

count <- dplyr::count
select <- dplyr::select

# Enable multi-threading for WGCNA to speed up network construction
allowWGCNAThreads(parallelly::availableCores())

# Parse the command-line argument to run the script for a specific stimulation
stim_i <- commandArgs(TRUE)[1]

# ------------------------------------------------------------------------------
# 2. Data Import & Filtering
# ------------------------------------------------------------------------------
txi <- read_rds("../1_quantification/results/gene_expression_txi.rds")

sample_table <- 
    tibble(sample_id = colnames(txi$counts)) |>
    extract(sample_id, c("group"), "[^_]+_(.+)", remove = FALSE) |>
    column_to_rownames("sample_id")

edger_res <- read_tsv("../2_dge/results/edger/results.tsv")
gene_names <- distinct(edger_res, gene_id, gene_name)  

# Restrict the network to genes identified as temporally significant 
# in the edgeR time-course analysis for this specific stimulation.
sig_genes <- 
    edger_res |>
    filter(stim == stim_i) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    pull(gene_id) |>
    unique()
   
# ------------------------------------------------------------------------------
# 3. Variance Stabilizing Transformation (VST)
# ------------------------------------------------------------------------------
# We use DESeq2's VST to stabilize variance prior to WGCNA.
dds <- 
    DESeqDataSetFromTximport(txi, sample_table, ~1) |>
    estimateSizeFactors() |>
    {function(x) x[rownames(x) %in% sig_genes, ]}() |>
    # Subset to only the Unstim baseline and the specific stimulation time points
    {function(x) x[, grepl(glue("_Unstim_0hrs|_{stim_i}_"), colnames(x))]}() |>
    # Require at least 4 counts in at least 3 samples
    {function(x) x[rowSums(counts(x, normalized = FALSE) >= 4 ) >= 3, ]}() |>
    vst()

# WGCNA requires genes as columns and samples as rows
count_matrix <- t(assay(dds))

# ------------------------------------------------------------------------------
# 4. Soft-Thresholding & Network Construction
# ------------------------------------------------------------------------------
# Explore multiple values of Beta (power) to find the one that best approximates 
# scale-free topology (where a few hub genes are highly connected, but most are not).

betas <- c(1:10, seq(12, 30, by = 2))

sft <- 
    pickSoftThreshold(count_matrix, 
		      powerVector = betas, 
		      networkType = "signed", 
		      blockSize = ncol(count_matrix),
		      verbose = 3)

## Plot beta vs r^2 between connectivity and frequency of connectivity
## and beta vs mean connectivity.
## The goal here is to find the minimum value of beta that maximizes r^2
## and minimizes connectivity.
#rsq_p <- 
#    ggplot(data = sft$fitIndices, 
#	   aes(x = Power, y = SFT.R.sq)) +
#    geom_hline(yintercept = 0.85, color = "red") +
#    geom_point() +
#    geom_line() +
#    scale_y_continuous(limits = c(0, 1),
#		       breaks = seq(0, 1, by = .2)) +
#    theme_bw()
#
#meank_p <-
#    ggplot(data = sft$fitIndices,
#	   aes(x = Power, y = mean.k.)) +
#    geom_point() +
#    geom_line() +
#    theme_bw()
#
#ggsave(glue("./plots/sft_{stim_i}.png"), 
#       rsq_p | meank_p, height = 3, width = 7)
#

merge_cut_height <- 0.2
detect_cut_height <- 0.99

# Default to a power of 30 if scale-free topology cannot be perfectly reached
beta_power <- ifelse(is.na(sft$powerEstimate), 30L, sft$powerEstimate)

# Temporarily reassign the cor function to WGCNA's optimized version to prevent bugs
cor <- WGCNA::cor

# Build the blockwise module network
network <- 
    blockwiseModules(count_matrix,
		     corType = "pearson",
		     networkType = "signed",
		     TOMType = "signed",
		     power = beta_power, 
		     maxBlockSize = ncol(count_matrix),
		     minModuleSize = 60,
		     mergeCutHeight = merge_cut_height,
		     detectCutHeight = detect_cut_height,
		     minKMEtoStay = 0.8,
		     reassignThreshold = 1e-3,
		     numericLabels = FALSE,
		     randomSeed = 1,
		     verbose = 3)

# Restore standard stats::cor function
cor <- stats::cor

module_genes <- enframe(network$colors, "gene_id", "module")

module_colors <- 
    module_genes |>
    count(module, sort = TRUE) |>
    filter(module != "grey") |>
    select(-n) |>
    deframe()

# ------------------------------------------------------------------------------
# 5. Network Visualizations (Dendrograms & Trees)
# ------------------------------------------------------------------------------
png(glue("./plots/dendro_{stim_i}.png"), 
    width = 10, height = 6, unit = "in", res = 200)
plotDendroAndColors(network$dendrograms[[1]], 
		    network$colors,
		    "merged",
		    dendroLabels = FALSE,
		    addGuide = TRUE,
		    hang = 0.03, 
		    guideHang = 0.05,
		    abHeight = detect_cut_height)
dev.off()

eigengene_tree <- 
    as.dist(1 - cor(select(network$MEs, -MEgrey))) |>
    hclust(method = "average")

png(glue("./plots/metree_{stim_i}.png"), 
    width = 8, height = 8, unit = "in", res = 300)
plot(eigengene_tree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = merge_cut_height, col = "red")
dev.off()

# ------------------------------------------------------------------------------
# 6. Intramodular Connectivity & kME calculations
# ------------------------------------------------------------------------------

# Calculate connectivity metrics to identify 'hub' genes driving each module
kim <-
    intramodularConnectivity.fromExpr(count_matrix, 
				      network$colors,
				      networkType = "signed",
				      power = beta_power,
				      scaleByMax = TRUE) 

kim_df <- 
    tibble(gene_id = colnames(count_matrix), kim = kim$kWithin) |>
    left_join(module_genes, join_by(gene_id)) |>
    left_join(gene_names) |>
    select(gene_id, gene_name, module, kim) |>
    filter(!is.na(kim))

# Calculate kME (Module Eigengene-based connectivity)
kme_all_df <- 
    signedKME(count_matrix, 
	      orderMEs(network$MEs), 
	      outputColumnName = "") |>
    as_tibble(rownames = "gene_id") |>
    select(-grey) |>
    pivot_longer(-gene_id, names_to = "module", values_to = "kme")

module_eigen <- as_tibble(network$MEs, rownames = "sample_name")

# ------------------------------------------------------------------------------
# 7. Trend Extraction & Plotting
# ------------------------------------------------------------------------------

# We select the top 500 genes based on kME rather than relying purely on 
# initial cluster assignment.
scaled_mean_counts <- 
    count_matrix |>
    scale() |>
    as_tibble(rownames = "sample_name") |>
    separate(sample_name, c("subject_id", "stim", "timep"), sep = "_") |>
    mutate(timep = sub("hrs$", "", timep),
	   timep = factor(as.integer(timep)),
	   stim = factor(stim, levels = c("Unstim", stim_i))) |>
    pivot_longer(-c(subject_id, stim, timep), names_to = "gene_id") |>
    group_by(stim, timep, gene_id) |>
    summarise(mean_scaled_counts = mean(value)) |>
    ungroup()

# Extract top 500 genes per module based on eigengene connectivity
top_kme <- 
    kme_all_df |>
    filter(kme >= 0.9) |>
    group_by(gene_id) |>
    slice_max(kme) |>
    ungroup() |>
    group_by(module) |>
    top_n(500, kme) |>
    ungroup() |>
    arrange(module, desc(kme))

top_kme_counts <- 
    top_kme |>
    left_join(scaled_mean_counts, join_by(gene_id)) |>
    mutate(module = factor(module, levels = module_colors))

mod_summ <- 
    module_eigen |>
    column_to_rownames("sample_name") |>
    scale() |>
    as_tibble(rownames = "sample_name") |>
    pivot_longer(-sample_name, names_to = "module") |>
    mutate(module = str_remove(module, "^ME")) |>
    separate(sample_name, c("donor_id", "stim", "timep"), sep = "_") |>
    group_by(timep, module) |>
    summarise(value = mean(value)) |>
    ungroup() |>
    filter(module != "grey") |>
    mutate(timep = str_remove(timep, "hrs$"),
	   timep = factor(as.integer(timep)),
	   module = factor(module, levels = module_colors))

trends_plot <- 
    ggplot(data = top_kme_counts,
	   aes(x = timep, y = mean_scaled_counts)) +
    geom_line(aes(group = gene_id, color = module), 
	      alpha = .05) +
    geom_line(data = mod_summ, 
	      aes(x = timep, y = value, group = module)) +
    geom_vline(xintercept = unique(top_kme_counts$timep), 
	       linetype = 2, linewidth = .1) +
    scale_color_manual(values = setNames(module_colors, module_colors)) +
    facet_wrap(~module, scale = "free", nrow = 2) +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  legend.position = "none") +
    labs(x = "Hours", y = "Average scaled expression")

ggsave(glue("./plots/trends_{stim_i}.png"), trends_plot, 
	    height = 3.5, width = 7)
   
# ------------------------------------------------------------------------------
# 8. Gene Ontology (GO) Enrichment
# ------------------------------------------------------------------------------
run_enrichment <- function(gene_list, background_list) {

    enrichGO(gene = gene_list,
	     universe = background_list,
	     OrgDb = org.Hs.eg.db,
	     ont = "BP",
	     keyType = "ENSEMBL",
	     pAdjustMethod = "fdr",
	     qvalueCutoff = 0.1,
	     readable = TRUE)
}

# The background universe consists of all genes that entered WGCNA
bkg_genes <- 
    kme_all_df |>
    distinct(gene_id) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    pull(gene_id)

plan(multisession, workers = availableCores())

go_res <- 
    top_kme |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    {function(x) split(x, x$module)}() |>
    map("gene_id") |>
    future_map(run_enrichment, background_list = bkg_genes)

plan(sequential)

go_res_df <- go_res |>
    map_df(as_tibble, .id = "module") 

go_res_top <- go_res_df |>
    mutate(Description = str_trunc(Description, width = 36)) |>
    group_by(module) |>
    top_n(10, -log10(pvalue)) |>
    ungroup()

go_plot <- 
    ggplot(data = go_res_top, 
       aes(x = "1", 
	   y = reorder_within(Description, by = -log10(pvalue), within = module))) +
    geom_point(aes(fill = -log10(pvalue), size = Count), shape = 21, stroke = .25) +
    scale_y_reordered() +
    scale_fill_gradient(low = "beige", high = "tomato4") +
    facet_wrap(~factor(module, levels = module_colors),
	       nrow = 2, scale = "free_y", drop = FALSE) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
	  axis.text.y = element_text(size = 9, margin = margin(r = -2, unit = "lines")),
	  axis.title = element_blank(),
	  panel.grid = element_blank(),
	  legend.position = "top",
	  legend.title = element_text(size = 10),
	  legend.text = element_text(size = 10),
	  plot.background = element_rect(fill = "white", color = "white"),
	  strip.text = element_text(size = 10, hjust = .5),
	  strip.clip = "off") +
    labs(fill = expression("-log"["10"]("p"))) +
    guides(fill = guide_colorbar(barheight = .7, barwidth = 10),
	   size = guide_legend(override.aes = list(fill = "black")))

ggsave(glue("./plots/go_{stim_i}.png"), 
       go_plot, width = 11, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# 9. Data Export
# ------------------------------------------------------------------------------
signedKME(count_matrix, orderMEs(network$MEs), outputColumnName = "") |>
    as_tibble(rownames = "gene_id") |>
    write_tsv(glue("data/{stim_i}_kme.tsv"))

write_rds(network, glue("data/{stim_i}_network.rds"))
write_tsv(module_genes, glue("./data/{stim_i}_modules.tsv"))
write_tsv(module_eigen, glue("./data/{stim_i}_eigen.tsv"))
write_rds(count_matrix, glue("./data/{stim_i}_counts.rds"))
write_tsv(kim_df, glue("./data/{stim_i}_kim.tsv"))
write_tsv(go_res_df, glue("./data/{stim_i}_go.tsv"))
