# Data handling
library(tidyverse)
library(glue)
library(furrr)

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

# Set WGCNA cores
allowWGCNAThreads(parallelly::availableCores())

# Expression data
dat <- read_rds("./data/gene_expression.rds")

sample_table <- 
    tibble(sample_id = colnames(dat$counts)) |>
    extract(sample_id, c("group"), "[^_]+_(.+)", remove = FALSE) |>
    column_to_rownames("sample_id")

edger_res <- read_tsv("../results/edger/results.tsv")
gene_names <- distinct(edger_res, gene_id, gene_name)  

## Use genes that are significant in the time course
stim_i <- commandArgs(TRUE)[1]

sig_genes <- 
    edger_res |>
    filter(stim == stim_i) |>
    pull(gene_id) |>
    unique()
    
## Create DESeq2 object for normalization of expression values.
## ~1 indicates no design, since we will not perform differential gene expression analysis here.
## Normalize counts using the variance-stabilizing transformation (VST)
dds <- 
    DESeqDataSetFromTximport(dat, sample_table, ~1) |>
    estimateSizeFactors() |>
    {function(x) x[rownames(x) %in% sig_genes, ]}() |>
    {function(x) x[, grepl(glue("_Unstim_0hrs|_{stim_i}_"), colnames(x))]}() |>
    {function(x) x[rowSums(counts(x, normalized = FALSE) >= 4 ) >= 3, ]}() |>
    vst()

## Transpose expression matrix to use with WGCNA
count_matrix <- t(assay(dds))

# Run WGCNA

## Explore multiple values of Beta (power) to find the one that satisfies scale-free topology
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
## That's because we believe that, in gene networks, there should be few genes
## with many connections (hub genes), and most genes should connect to a small
## number of genes.
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
### Build WGCNA network 
### we need to reassing the 'cor()' function to avoid a bug in WGCNA
merge_cut_height <- 0.2
detect_cut_height <- 0.99

beta_power <- ifelse(is.na(sft$powerEstimate), 30L, sft$powerEstimate)

cor <- WGCNA::cor

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

cor <- stats::cor

module_genes <- enframe(network$colors, "gene_id", "module")

module_colors <- 
    module_genes |>
    count(module, sort = TRUE) |>
    filter(module != "grey") |>
    select(-n) |>
    deframe()

## plot dendogram to visualize modules
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

# plot eigengene tree
eigengene_tree <- 
    as.dist(1 - cor(select(network$MEs, -MEgrey))) |>
    hclust(method = "average")

png(glue("./plots/metree_{stim_i}.png"), 
    width = 8, height = 8, unit = "in", res = 300)
plot(eigengene_tree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = merge_cut_height, col = "red")
dev.off()


# Gene classification into modules
## Intramodule connectivity
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

## kME
kme_all_df <- 
    signedKME(count_matrix, 
	      orderMEs(network$MEs), 
	      outputColumnName = "") |>
    as_tibble(rownames = "gene_id") |>
    select(-grey) |>
    pivot_longer(-gene_id, names_to = "module", values_to = "kme")


## Module eigengenes
## A module eigengene is like a hypothetical gene the represents the module
module_eigen <- as_tibble(network$MEs, rownames = "sample_name")

## Plot trends
### Module assignment and higher kME are not always concordant.
### While module eigengene is based on original assigments (based on adjacency not correlation).
### Some genes have very low kME with the assigned module, and the threshold used in blockwideModules is not working properly.
### I'm selecting top 500 genes by highest kME with a module to test with GO
### By suggestion of WGCNA main author.
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

kme_df <- 
    kme_all_df |>
    filter(kme >= 0.9) |>
    group_by(gene_id) |>
    slice_max(kme) |>
    ungroup()

top_kme <- kme_df |>
    group_by(module) |>
    top_n(500, kme) |>
    ungroup() |>
    arrange(module, desc(kme))

top_kme_counts <- 
    top_kme |>
    left_join(scaled_mean_counts, join_by(gene_id)) |>
    mutate(module = factor(module, levels = module_colors))

mod_summ <- module_eigen |>
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
   
# GO
#library(fgsea)
#
#pathwaysgo <- gmtPathways("/lab-share/IM-Gutierrez-e2/Public/References/msigdb/c5.all.v2022.1.Hs.symbols.gmt.txt")
#gobp <- keep(pathwaysgo, grepl("^GOBP", names(pathwaysgo)))
#
#filtered_genes <- 
#    gene_names |>
#    add_count(gene_name) |>
#    filter(n == 1) |>
#    select(gene_id, gene_name)
#
#gene_list <-
#    kme_all_df |>
#    inner_join(filtered_genes) |>
#    select(module, gene_name, kme) |>
#    {function(x) split(x, x$module)}() |>
#    map(~select(., -module)) |>
#    map(~arrange(., desc(kme))) |>
#    map(deframe)
#
#gsea_res <- 
#    map(gene_list, 
#	~fgsea(pathways = gobp, stats = ., nproc = future::availableCores()))
#
#gsea_df <- 
#    bind_rows(gsea_res, .id = "module") |>
#    as_tibble()
#

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


## Run GO
## it takes a few minutes
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


# Plot correlation between Lupus-associated genes and modules in a heatmap
# List of Lupus-associated genes
sle_genes <- 
    c("PTPN22", 
      "FCGR2A", 
      "TNFSF4", 
      "NCF2", 
      "IL10",
      "LYST",
      "SPRED2", 
      "IFIH1", 
      "STAT1", "STAT4", 
      "IKZF2", 
      "PXK", 
      "IL12A", 
      "BANK1", 
      "TCF7", "SKP1",
      "TNIP1", 
      "MIR3142HG", 
      "UHRF1BP1", 
      "ATG5", "PRDM1",
      "TNFAIP3", 
      "JAZF1",
      "IKZF1", 
      "IRF5", "TNPO3",
      "BLK", 
      "WDFY4", 
      "ARID5B", 
      "IRF7", 
      "CD44",
      "DHCR7", "NADSYN1",
      "ETS1", "FLI1",
      "SH2B3", 
      "SLC15A4",
      "RAD51B",
      "CSK", 
      "CIITA", "SOCS1", 
      "ITGAM", "ITGAX", 
      "IRF8",
      "PLD2",
      "IKZF3", 
      "TYK2", 
      "UBE2L3", 
      "TASL",
      "IRAK1", "MECP2",
      "IKBKE", "IL10" 
    )


sle_genes_cormatrix <- 
    gene_names |>
    filter(gene_name %in% sle_genes) |>
    left_join(kme_all_df, join_by(gene_id)) |>
    filter(!is.na(module)) |>
    #group_by(gene_id) |>
    #filter(any(kme >= 0.85)) |>
    #ungroup() |>
    select(-gene_id) |>
    pivot_wider(names_from = module, values_from = kme) |>
    column_to_rownames("gene_name") |>
    select(all_of(module_colors)) |>
    {function(x) setNames(x, glue("M{1:ncol(x)}"))}() |>
    data.matrix()

png(glue("./plots/slegenes_{stim_i}.png"), 
    units = "in", height = 6, width = 4, res = 400)
pheatmap(sle_genes_cormatrix,
	 fontsize = 9, angle_col = 0, cluster_rows = TRUE, cluster_cols = FALSE, 
	 color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100),
	 display_numbers = round(sle_genes_cormatrix, 2),
	 legend = FALSE)
dev.off()

# Save data for analysis of module preservation across stims
signedKME(count_matrix, orderMEs(network$MEs), outputColumnName = "") |>
    as_tibble(rownames = "gene_id") |>
    write_tsv(glue("data/{stim_i}_kme.tsv"))

write_rds(network, glue("data/{stim_i}_network.rds"))
write_tsv(module_genes, glue("./data/{stim_i}_modules.tsv"))
write_tsv(module_eigen, glue("./data/{stim_i}_eigen.tsv"))
write_rds(count_matrix, glue("./data/{stim_i}_counts.rds"))
write_tsv(kim_df, glue("./data/{stim_i}_kim.tsv"))
write_tsv(go_res_df, glue("./data/{stim_i}_go.tsv"))
#
#
### TOM
##tom <- 
##    TOMsimilarityFromExpr(count_matrix, 
##			  corType = "pearson",
##			  networkType = "signed",
##			  power = beta_power,
##			  nThreads = availableCores())
##
##colnames(tom) <- rownames(tom) <- colnames(count_matrix)
##
##top_10_hub <- 
##    kim_df |>
##    group_by(module) |>
##    top_n(10, kim) |>
##    ungroup()
##
##tom_hub <- 
##    tom[top_10_hub$gene_id, top_10_hub$gene_id] |>
##    as_tibble(rownames = "from_id") |>
##    pivot_longer(-from_id, names_to = "to_id") |>
##    left_join(select(top_10_hub, -kim), join_by(from_id == "gene_id")) |>
##    left_join(select(top_10_hub, gene_id, gene_name), join_by(to_id == "gene_id")) |>
##    select(from_id = gene_name.x, to_id = gene_name.y, module, value) |>
##    filter(from_id != to_id) |>
##    arrange(desc(value))
##
##write_tsv(tom_hub, glue("./data/{stim_i}_TOM_hub.tsv"))
##
