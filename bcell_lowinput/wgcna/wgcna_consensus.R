# Data handling
library(tidyverse)
library(glue)

# import data
library(tximport)
library(DESeq2)

# Run main analyses
library(WGCNA)

# Gene Ontology analysis
library(clusterProfiler)
library(org.Hs.eg.db)

# Plotting
library(pheatmap)
library(RColorBrewer)
library(tidytext)
library(patchwork)

count <- dplyr::count
select <- dplyr::select

# Expression data
dat <- read_rds("./data/gene_expression.rds")

sample_table <- 
    tibble(sample_id = colnames(dat$counts)) |>
    extract(sample_id, c("group"), "[^_]+_(.+)", remove = FALSE) |>
    column_to_rownames("sample_id")

edger_res <- read_tsv("../results/edger/results.tsv")

# Run WGCNA
allowWGCNAThreads()

stim_1 <- "BCR"
stim_2 <- "TLR7"
stim_3 <- "BCR-TLR7"

## Use genes that are significant in the time course
sig_genes_1 <- edger_res |>
    filter(stim == stim_1, FDR < 0.05) |>
    pull(gene_id)

sig_genes_2 <- edger_res |>
    filter(stim == stim_2, FDR < 0.05) |>
    pull(gene_id)

sig_genes_3 <- edger_res |>
    filter(stim == stim_3, FDR < 0.05) |>
    pull(gene_id)
    
sig_genes <- unique(c(sig_genes_1, sig_genes_2, sig_genes_3))

## Create DESeq2 object for normalization of expression values.
## ~1 indicates no design, since we will not perform differential gene expression analysis here.
## Select samples with library size > 2MM reads.
## Select the same genes for all conditions (requirement of consensus analysis).
## Filter for genes with normalized counts > 10 in at least 3 samples.
## Normalize counts using the variance-stabilizing transformation (VST)

dds <- 
    DESeqDataSetFromTximport(dat, sample_table, ~1) |>
    estimateSizeFactors() |>
    {function(x) x[, colSums(counts(x)) > 2e6]}() |>
    {function(x) x[rownames(x) %in% sig_genes, ]}() |>
    {function(x) x[, grepl(glue("_Unstim_0hrs|_{stim_1}_|_{stim_2}_|_{stim_3}_"), colnames(x))]}() |>
    {function(x) x[rowSums(counts(x, normalized = TRUE) >= 10 ) >= 3, ]}() |>
    vst()

count_matrix_1 <- dds |> 
    {function(x) x[, grepl(glue("_Unstim_0hrs|_{stim_1}_"), colnames(x))]}() |>
    assay() |>
    t()

count_matrix_2 <- dds |> 
    {function(x) x[, grepl(glue("_Unstim_0hrs|_{stim_2}_"), colnames(x))]}() |>
    assay() |>
    t()

count_matrix_3 <- dds |> 
    {function(x) x[, grepl(glue("_Unstim_0hrs|_{stim_3}_"), colnames(x))]}() |>
    assay() |>
    t()

good_genes_1 <- which(apply(count_matrix_1, 2, var) > 0)
good_genes_2 <- which(apply(count_matrix_2, 2, var) > 0)
good_genes_3 <- which(apply(count_matrix_3, 2, var) > 0)

good_genes <- 
    intersect(good_genes_1, good_genes_2) |>
    intersect(good_genes_3)

count_matrix_1 <- count_matrix_1[, good_genes]
count_matrix_2 <- count_matrix_2[, good_genes]
count_matrix_3 <- count_matrix_3[, good_genes]


## Explore multiple values of Beta (power) to find the one that satisfies scale-free topology
betas <- c(1:10, seq(12, 30, by = 2))

sft_1 <- 
    pickSoftThreshold(count_matrix_1, 
		      powerVector = betas, 
		      networkType = "signed", 
		      blockSize = ncol(count_matrix_1),
		      verbose = 3)

sft_2 <- 
    pickSoftThreshold(count_matrix_2, 
		      powerVector = betas, 
		      networkType = "signed", 
		      blockSize = ncol(count_matrix_2),
		      verbose = 3)

sft_3 <- 
    pickSoftThreshold(count_matrix_3, 
		      powerVector = betas, 
		      networkType = "signed", 
		      blockSize = ncol(count_matrix_3),
		      verbose = 3)

best_betas <- 
    list(sft_1, sft_2, sft_3) |>
    map_int("powerEstimate")


## Build WGCNA network 
## we need to reassing the 'cor()' function to avoid a bug in WGCNA
cut_height <- 0.3

multi_expr <- list()
multi_expr[[1]] <- list(data = count_matrix_1)
multi_expr[[2]] <- list(data = count_matrix_2)
multi_expr[[3]] <- list(data = count_matrix_3)
names(multi_expr) <- c(stim_1, stim_2, stim_3)

cor <- WGCNA::cor

network <- 
    blockwiseConsensusModules(multi_expr,
			      corType = "bicor",
			      networkType = "signed",
			      TOMType = "signed",
			      power = best_betas, 
			      maxBlockSize = 30000,
			      minModuleSize = 30, 
			      mergeCutHeight = cut_height,
			      detectCutHeight = 0.995,
			      randomSeed = 1,
			      verbose = 3)

cor <- stats::cor


png(glue("./plots/wgcna_cons_{stim_1}_{stim_2}_{stim_3}_dendro.png"), 
    width = 10, height = 6, unit = "in", res = 200)
plotDendroAndColors(network$dendrograms[[1]], 
		    network$colors,
		    "merged",
		    dendroLabels = FALSE,
		    addGuide = TRUE,
		    hang = 0.03, 
		    guideHang = 0.05)
dev.off()

## Gene classification into modules
module_genes <- enframe(network$colors, "gene_id", "module")

## Module eigengenes
## A module eigengene is like a hypothetical gene the represents the module
module_eigen <- 
    network$multiMEs |> 
    map("data") |>
    map_df(~as_tibble(., rownames = "sample_name"), .id = "stim")

# plot trends
count_list <-
    list(count_matrix_1, count_matrix_2, count_matrix_3) |>
    setNames(c(stim_1, stim_2, stim_3))

scaled_mean_counts <- count_list |>
    map_df(~scale(.) |>
	   as_tibble(rownames = "sample_name") |>
	   separate(sample_name, c("subject_id", "stim", "timep"), sep = "_") |>
	   mutate(timep = sub("hrs$", "", timep),
		  timep = factor(as.integer(timep))) |>
	   pivot_longer(-c(subject_id, stim, timep), names_to = "gene_id") |>
	   group_by(stim, timep, gene_id) |>
	   summarise(mean_scaled_counts = mean(value)) |>
	   ungroup(), .id = "stim")

kme_df <-
    map_df(setNames(1:3, names(count_list)),
	   ~signedKME(count_list[[.x]], orderMEs(network$multiMEs[[.x]]$data)) |>
	   as_tibble(rownames = "gene_id") |>
	   select(-kMEgrey) |>
	   pivot_longer(-gene_id, names_to = "module"), 
	   .id = "stim")

# get average module membership for genes across stims
top_kme <- kme_df |>
    group_by(module, gene_id) |>
    summarise(value = mean(value)) |>
    group_by(gene_id) |>
    slice_max(value) |>
    group_by(module) |>
    top_n(500, abs(value)) |>
    ungroup() |>
    arrange(module, desc(abs(value)))

mod_summ <- module_eigen |>
    {function(x) split(x, x$stim)}() |>
    map_df(~select(., -stim) |>
	   column_to_rownames("sample_name") |>
	   scale() |>
	   as_tibble(rownames = "sample_name") |>
	   pivot_longer(-sample_name, names_to = "module") |>
	   mutate(module = str_remove(module, "^ME")) |>
	   separate(sample_name, c("donor_id", "stim", "timep"), sep = "_") |>
	   group_by(timep, module) |>
	   summarise(value = mean(value)) |>
	   ungroup() |>
	   mutate(timep = str_remove(timep, "hrs$"),
		  timep = factor(as.integer(timep))) |>
	   filter(module != "grey"),
	   .id = "stim")

module_colors <- top_kme |> 
    distinct(module) |>
    mutate(module = str_remove(module, "kME")) |>
    pull(module) |>
    {function(x) setNames(x, x)}()

module_colors["yellow"] <- "goldenrod3"

trends_plot <- 
    ggplot() + 
    geom_line(data = mod_summ, 
	      aes(x = timep, y = value, 
		  linetype = stim, color = module,
		  group = interaction(module, stim))) +
    geom_vline(xintercept = unique(mod_summ$timep), 
	       linetype = 2, linewidth = .1) +
    scale_color_manual(values = module_colors, guide = NULL) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    facet_wrap(~module, scale = "free") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Hours", y = "Average scaled expression")

ggsave(glue("./plots/wgcna_cons_{stim_1}_{stim_2}_{stim_3}_trends.png"), trends_plot, 
	    height = 3.5, width = 7)

