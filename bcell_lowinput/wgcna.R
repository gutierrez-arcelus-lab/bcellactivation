library(tidyverse)
library(tximport)
library(DESeq2)
library(WGCNA)
library(patchwork)

# Sample meta data
meta_data <- read_tsv("./data/sample_decode.tsv") |>
    separate(sample_name, c("name", "stim", "time"), sep = "_", remove = FALSE) |>
    mutate(time = factor(time, levels = str_sort(unique(time), numeric = TRUE))) |>
    arrange(stim, time, name)

sample_table <- meta_data |>
    unite("group", c(stim, time), sep = ".") |>
    select(sample_id, group) |>
    column_to_rownames("sample_id")

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
    sprintf("./results/salmon/%s/quant.sf", meta_data$sample_id) |>
    setNames(meta_data$sample_id)

txi <- tximport(salmon_files, 
		type = "salmon", 
		tx2gene = select(tx_to_gene, tx_id, gene_id))

dds <- DESeqDataSetFromTximport(txi, sample_table, ~1)
dds <- estimateSizeFactors(dds)
keep_y <- rowSums(counts(dds, normalized = TRUE) >= 10 ) >= 3
dds <- dds[keep_y, ]
dds <- collapseReplicates(dds, meta_data$sample_name, meta_data$sample_id)

gsg <- goodSamplesGenes(t(counts(dds)))
gsg$allOK

dds <- vst(dds)

# WGCNA
allowWGCNAThreads()

vst_counts <- t(assay(dds))

powers <- c(1:10, seq(12, 50, by = 2))

sft <- 
    pickSoftThreshold(vst_counts, 
		      powerVector = powers, 
		      networkType = "signed", 
		      blockSize = 30000,
		      verbose = 5)

rsq_p <- 
    ggplot(data = sft$fitIndices, 
	   aes(x = Power, y = SFT.R.sq)) +
    geom_point() +
    geom_line() +
    theme_bw()

meank_p <-
    ggplot(data = sft$fitIndices,
	   aes(x = Power, y = mean.k.)) +
    geom_point() +
    geom_line() +
    theme_bw()

ggsave("./plots/wgcna_powers.png", rsq_p | meank_p, height = 3)

cor <- WGCNA::cor

bwnet <- 
    blockwiseModules(vst_counts,
		     corType = "bicor",
		     networkType = "signed",
		     power = 14, 
		     maxBlockSize = 30000,
		     minModuleSize = 30,
		     mergeCutHeight = .25,
		     numericLabels = FALSE,
		     randomSeed = 1,
		     verbose = 3)

cor <- stats::cor

mod_eigengenes <- as_tibble(bwnet$MEs, rownames = "sample_name")

mod_gene <- enframe(bwnet$colors, "gene_id", "module")

mod_gene |> dplyr::count(module, sort = TRUE)

png("./plots/wgcna_dendrogram.png", width = 8, height = 4, unit = "in", res = 300)
plotDendroAndColors(bwnet$dendrograms[[1]], 
		    bwnet$colors,
		    "merged",
		    dendroLabels = FALSE,
		    addGuide = TRUE,
		    hang = 0.03, 
		    guideHang = 0.05)
dev.off()

mod_eigen_df <- mod_eigengenes |>
    pivot_longer(-sample_name, names_to = "module") |>
    separate(sample_name, c("subject_id", "stim", "timep"), sep = "_") |>
    mutate(timep = sub("hrs$", "", timep),
	   timep = factor(as.integer(timep)),
	   stim = factor(stim, levels = c("Unstim", "IL4", "CD40L", "BCR", "TLR7", "TLR9", "BCR-TLR7", "DN2"))) |>
    group_by(stim, timep, module) |>
    summarise(value = mean(value)) |>
    ungroup()

tp <- 
    ggplot(data = mod_eigen_df,
	   aes(x = timep, y = module)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(high = "tomato3", mid = "white", low = "midnightblue",
			 guide = guide_colorbar(barheight = .5, barwidth = 10)) +
    facet_grid(~stim, scale = "free_x", space = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Hours", fill = NULL)

ggsave("./plots/wgcna_tp.png", tp, height = 3, width = 10)

count_norm_df <- vst_counts |>
    as_tibble(rownames = "sample_name") |>
    separate(sample_name, c("subject_id", "stim", "timep"), sep = "_") |>
    mutate(timep = sub("hrs$", "", timep),
	   timep = factor(as.integer(timep)),
	   stim = factor(stim, levels = c("Unstim", "IL4", "CD40L", "BCR", "TLR7", "TLR9", "BCR-TLR7", "DN2"))) |>
    pivot_longer(-c(subject_id, stim, timep), names_to = "gene_id") |>
    left_join(mod_gene, join_by(gene_id)) |>
    left_join(distinct(tx_to_gene, gene_id, gene_name), join_by(gene_id)) |>
    select(subject_id, stim, timep, gene_id, gene_name, module, value)

count_mod_summ <- count_norm_df |>
    group_by(stim, timep, gene_id, gene_name, module) |>
    summarise(value = mean(value)) |>
    ungroup()

yellow_p <- 
    ggplot(data = count_mod_summ |> filter(module == "yellow"),
       aes(x = timep, y = value)) +
    geom_line(aes(group = gene_id), alpha = .25) +
    facet_grid(~stim, scale = "free_x", space = "free") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Hours")

ggsave("./plots/wgcna_yellow.png", yellow_p, height = 3, width = 12)

black_p <- 
    ggplot(data = count_mod_summ |> filter(module == "black"),
       aes(x = timep, y = value)) +
    geom_line(aes(group = gene_id), alpha = .1) +
    facet_grid(~stim, scale = "free_x", space = "free") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Hours")

ggsave("./plots/wgcna_black.png", yellow_p, height = 3, width = 12)


count_mod_summ |> filter(gene_name == "TNFSF4") |> print(n = Inf)

# VST transform seems to bring very different expression levels close
# so samples don't look different anymore.
# Try TPM?
count_mod_summ |> filter(gene_name == "STAT1") |> print(n = Inf)


# Try TPM
keep_samples <- txi$counts |> 
    colSums() |>
    enframe("sample_id", "lib_size") |>
    left_join(meta_data, join_by(sample_id)) |>
    group_by(sample_name) |>
    slice_max(lib_size) |>
    ungroup()

tpm_matrix <- txi$counts[, keep_samples$sample_id]
colnames(tpm_matrix) <- keep_samples$sample_name
keep_tpm <- rowSums(tpm_matrix >= 10 ) >= 3
tpm_matrix <- tpm_matrix[keep_tpm, ]
tpm_matrix <- t(tpm_matrix)

powers <- c(1:10, seq(12, 50, by = 2))

sft <- 
    pickSoftThreshold(tpm_matrix, 
		      powerVector = powers, 
		      networkType = "signed", 
		      blockSize = 30000,
		      verbose = 5)

rsq_p <- 
    ggplot(data = sft$fitIndices, 
	   aes(x = Power, y = SFT.R.sq)) +
    geom_point() +
    geom_line() +
    theme_bw()

meank_p <-
    ggplot(data = sft$fitIndices,
	   aes(x = Power, y = mean.k.)) +
    geom_point() +
    geom_line() +
    theme_bw()

ggsave("./plots/wgcna_powers_tpm.png", rsq_p | meank_p, height = 3)

cor <- WGCNA::cor

bwnet_tpm <- 
    blockwiseModules(tpm_matrix,
		     corType = "bicor",
		     networkType = "signed",
		     power = 16, 
		     maxBlockSize = 30000,
		     minModuleSize = 30,
		     mergeCutHeight = .25,
		     numericLabels = FALSE,
		     randomSeed = 1,
		     verbose = 3)

cor <- stats::cor

mod_eigengenes_tpm <- as_tibble(bwnet_tpm$MEs, rownames = "sample_name")

mod_gene_tpm <- enframe(bwnet_tpm$colors, "gene_id", "module")

mod_gene_tpm |> dplyr::count(module, sort = TRUE)

png("./plots/wgcna_dendrogram_tpm.png", width = 8, height = 4, unit = "in", res = 300)
plotDendroAndColors(bwnet_tpm$dendrograms[[1]], 
		    bwnet_tpm$colors,
		    "merged",
		    dendroLabels = FALSE,
		    addGuide = TRUE,
		    hang = 0.03, 
		    guideHang = 0.05)
dev.off()

mod_eigen_tpm_df <- mod_eigengenes_tpm |>
    pivot_longer(-sample_name, names_to = "module") |>
    separate(sample_name, c("subject_id", "stim", "timep"), sep = "_") |>
    mutate(timep = sub("hrs$", "", timep),
	   timep = factor(as.integer(timep)),
	   stim = factor(stim, levels = c("Unstim", "IL4", "CD40L", "BCR", "TLR7", "TLR9", "BCR-TLR7", "DN2"))) |>
    group_by(stim, timep, module) |>
    summarise(value = mean(value)) |>
    ungroup()

tp_tpm <- 
    ggplot(data = mod_eigen_tpm_df,
	   aes(x = timep, y = module)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(high = "tomato3", mid = "white", low = "midnightblue",
			 guide = guide_colorbar(barheight = .5, barwidth = 10)) +
    facet_grid(~stim, scale = "free_x", space = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Hours", fill = NULL)

ggsave("./plots/wgcna_tp_tpm.png", tp, height = 3, width = 10)


