library(DESeq2)
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(tibble)
library(forcats)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(cowplot)

stim_order <- 
    c(sprintf("unst %s", c(0, 24)), 
      sprintf("IL4 %s", c(24, 72)),
      sprintf("BCR %s", c(24, 72)),
      sprintf("TLR7 %s", c(24, 72)),
      sprintf("DN2 %s", c(24, 72))
    )
      
stim_colors <- 
    c("grey95", "grey70",
      "goldenrod1", "goldenrod4",
      brewer.pal(n = 9, "Blues")[c(3, 8)],
      brewer.pal(n = 9, "Greens")[c(3, 8)],
      paste0("tomato", c(2, 4))
      )

names(stim_colors) <- stim_order

donor_ids <- read_csv("./samplesheet.csv") |>
    mutate(donor_id = basename(fastq_1)) |>
    extract(donor_id, "donor_id", "[^_]+_([^_]+)_.+") |>
    mutate(replicate = paste0("REP", replicate),
	   stim = sub("_", " ", sample)) |>
    select(stim, donor_id, donor = replicate) |>
    distinct()

pca_df <- 
    "./results/bwa/merged_library/macs2/narrow_peak/consensus/deseq2/consensus_peaks.mLb.clN.pca.vals.txt" |>
    read_tsv() |>
    extract(sample, c("stim", "donor"), "(.+_.+)_(.+)") |>
    mutate(stim = sub("_", " ", stim)) |>
    left_join(donor_ids) |>
    mutate(stim = factor(stim, levels = stim_order))

pca_plot <- ggplot(pca_df, aes(`PC1: 26% variance`, `PC2: 15% variance`)) +
    geom_point(aes(fill = stim), shape = 21, size = 4) +
    geom_text_repel(aes(label = donor_id), size = 2.5, color = "grey30", min.segment.length = 0) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = "white"),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank())


ggsave("./plots/pca.png", pca_plot, width = 6, height = 5)

total_reads <- 
    "./results/multiqc/narrow_peak/B_cell_atac_nfcore_multiqc_report_data/multiqc_samtools_flagstat.txt" |>
    read_tsv() |>
    select(Sample, total_passed) |>
    extract(Sample, c("stim", "donor", "rep"), "([^_]+_\\d+)_(REP\\d)_(T\\d)") |>
    group_by(stim, donor) |>
    summarise(total_passed = sum(total_passed)) |>
    ungroup() |>
    mutate(stim = sub("_", " ", stim)) |>
    mutate(stim = factor(stim, levels = stim_order))

pca_plot_depth <- pca_df |>
    left_join(total_reads) |>
    ggplot(aes(`PC1: 26% variance`, `PC2: 15% variance`)) +
    geom_point(aes(fill = total_passed), shape = 21, size = 4) +
    geom_text_repel(aes(label = donor_id), size = 2.5, color = "grey30", min.segment.length = 0) +
    scale_fill_continuous(low = "lightyellow", high = "tomato3",
			  labels = scales::comma,
			  guide = guide_colorbar(barwidth = .5)) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = "white"),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank()) +
    labs(fill = "Total\nreads")

ggsave("./plots/pca_depth.png", pca_plot_depth, width = 6, height = 5)


mito <- 
    "./results/multiqc/narrow_peak/B_cell_atac_nfcore_multiqc_report_data/mqc_samtools-idxstats-mapped-reads-plot_Raw_Counts.txt" |>
    read_tsv() |>
    pivot_longer(-Sample, names_to = "chr", values_to = "n") |>
    extract(Sample, c("stim", "donor", "rep"), "([^_]+_\\d+)_(REP\\d)_(T\\d)") |>
    group_by(stim, donor, chr) |>
    summarise(n = sum(n)) |>
    group_by(stim, donor) |>
    mutate(p = n/sum(n)) |>
    ungroup() |>
    filter(chr == "chrM") |>
    select(-chr, -n) |>
    mutate(stim = sub("_", " ", stim)) |>
    mutate(stim = factor(stim, levels = stim_order))

pca_plot_mito <- pca_df |>
    left_join(mito) |>
    ggplot(aes(`PC1: 26% variance`, `PC2: 15% variance`)) +
    geom_point(aes(fill = p), shape = 21, size = 4) +
    geom_text_repel(aes(label = donor_id), size = 2.5, color = "grey30", min.segment.length = 0) +
    scale_fill_continuous(low = "white", high = "midnightblue",
			  labels = scales::comma,
			  guide = guide_colorbar(barwidth = .5)) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = "white"),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank()) +
    labs(fill = "% Mito")

ggsave("./plots/pca_mito.png", pca_plot_mito, width = 6, height = 5)


# try more peaks (5000 instead of 500)
pca_5000_df <- 
    "./results_deseq2/consensus_peaks.mLb.clN_pcadata_5000peaks.rds" |>
    read_rds() |>
    mutate(condition = sub("_", " ", condition)) |>
    extract(name, "donor", ".+(REP\\d)") |>
    select(stim = condition, donor, PC1, PC2) |>
    left_join(donor_ids) |>
    mutate(stim = factor(stim, levels = stim_order))

perc_var <-
    "./results_deseq2/consensus_peaks.mLb.clN_pcadata_5000peaks.rds" |>
    read_rds() |>
    attr("percentVar") |>
    {function(x) round(x * 100)}()

pca_5000_plot <- 
    ggplot(pca_5000_df, aes(PC1, PC2)) +
    geom_point(aes(fill = stim), shape = 21, size = 4) +
    geom_text_repel(aes(label = donor_id), size = 2.5, color = "grey30", min.segment.length = 0) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = "white"),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank()) +
    labs(x = sprintf("PC1: %s%% variance", perc_var[1]),
	 y = sprintf("PC2: %s%% variance", perc_var[2]))

ggsave("./plots/pca_5000VarPeaks.png", pca_5000_plot, width = 6, height = 5)



# DESeq2 PCA
library(DESeq2)
load("./results_deseq2/consensus_peaks.mLb.clN.dds.rld.RData")

# calculate the variance for each gene
rv <- rowVars(assay(rld))

# select the ntop genes by variance
select_top <- order(rv, decreasing = TRUE)[1:5000]

# perform a PCA on the data in assay(x) for the selected genes
pca_df <- as.data.frame(t(assay(rld)[select_top, ]))
pca <- prcomp(pca_df)

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

pca_loadings <- pca$rotation |> 
    as_tibble(rownames = "peak_id") |>
    select(peak_id, PC1, PC2) |>
    pivot_longer(-peak_id, names_to = "pc") |>
    group_by(pc, sign = value > 0) |>
    top_n(20, abs(value)) |>
    ungroup() |>
    left_join(select(peak_annot, peak_id, gene = `Gene Name`), by = "peak_id") |>
    mutate(peak_id = sub("Interval_", "", peak_id)) |>
    unite("peak_id", c(gene, peak_id))

pca_plot_2 <- ggplot(pca_loadings, 
		     aes(x = value, 
			 y = tidytext::reorder_within(peak_id, value, pc))) +
    geom_col(fill = "midnightblue", alpha = .5) +
    tidytext::scale_y_reordered() +
    facet_wrap(~pc, scales = "free_y") +
    labs(x = "Loading", y = NULL) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8))

ggsave("./plots/pca_loadings.png", pca_plot_2, width = 10)













## Annot 

peak_annot <- 
    "./results/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.annotatePeaks.txt" |> 
    read_tsv()

names(peak_annot)[1] <- "peak_id"

# DA

gwas_lead <- read_tsv("../colocalization/data/bentham_top_hits.tsv") |>
    mutate(region = as.character(region))

gwas_lead |> print(n = Inf)

coloc_res <- "../colocalization/results/bentham_region%d.tsv" |>
    sprintf(1:33) |>
    setNames(1:33) |>
    map_dfr(read_tsv, .id = "region") |>
    filter(h4 >= .90) |>
    left_join(gwas_lead, by = "region", multiple = "all") |>
    select(region, locus, study, gene_id, gene_name, molecular_trait_id, h4)


da_files <- 
    list.files("./results_deseq2", 
	       pattern = "\\.deseq2\\.FDR0\\.05\\.results\\.txt",
	       full.names = TRUE)

da_names <- 
    basename(da_files) |>
    strsplit("\\.") |>
    map_chr(1)

da_all <- da_files |>
    setNames(da_names) |>
    {function(x) x[!grepl("72", x)]}() |>
    map_dfr(~read_tsv(.) |> select(Geneid:padj), .id = "contrast")

da_all_tmp <- da_all |>
    dplyr::count(contrast) |>
    separate(contrast, c("c0", "c1"), sep = "vs")

stims <- c("unst_0", "unst_24", "IL4_24", "TLR7_24", "BCR_24", "DN2_24")

da_all_summ <- 
    bind_rows(da_all_tmp, select(da_all_tmp, c0 = c1, c1 = c0, n)) |>
    complete(c0, c1, fill = list(n = 0)) |>
    mutate_at(vars(c0:c1), ~factor(., levels = stims)) |>
    mutate(c0l = as.integer(c0), 
	   c1l = as.integer(c1)) |>
    mutate(cmin = pmin(c0l, c1l), 
	   cmax = pmax(c0l, c1l),
	   c0 = stims[cmin],
	   c1 = stims[cmax]) |>
    distinct(cmin, cmax, .keep_all = TRUE) |>
    arrange(cmin, cmax) |>
    select(c0, c1, n) |>
    mutate_at(vars(c0:c1), fct_inorder) |>
    mutate(n = ifelse(n == 0 & c0 == c1, NA, n))

da_summ_plot <- 
    ggplot(da_all_summ, aes(c0, c1)) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = scales::comma(n)), 
	      color = "white", fontface = "bold") +
    scale_x_discrete(labels = function(x) sub("_", " ", x)) +
    scale_y_discrete(labels = function(x) sub("_", " ", x)) +
    scale_fill_gradient(low = "midnightblue", high = "tomato3", 
			na.value = "white") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = "white"),
	  plot.margin = margin(1, 1, 1, 1, unit = "cm"),
	  panel.grid = element_blank(),
	  axis.title = element_blank(),
	  axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 1),
	  axis.text.y = element_text(size = 12),
	  legend.position = "none")

ggsave("./plots/diffaccess_summ.png", da_summ_plot, 
       height = 5.5, width = 5.5, dpi = 600)

pca_out <- 
    plot_grid(pca_5000_plot, NULL, da_summ_plot,
	      nrow = 1, rel_widths = c(1, .1, 1),
	      labels = c("A)", "", "B)")) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/pca_da.png", width = 10, height = 5)





## stim vs IL4
da_bcr24_il24 <- "./results_deseq2/BCR_24vsIL4_24.deseq2.FDR0.01.results.txt" |>
    read_tsv() |>
    arrange(padj)

da_dn24_il24 <- "./results_deseq2/DN2_24vsIL4_24.deseq2.FDR0.01.results.txt" |>
    read_tsv() |>
    arrange(padj)

gwas_region <- gwas_lead |>
    mutate(left = pos - 2.5e5, right = pos + 2.5e5)

da_unst24_dn24 |>
    select(Geneid, Chr, Start, End, log2FoldChange, padj) |>
    filter(Geneid %in% da_dn24_il24$Geneid) |>
    inner_join(gwas_region, join_by(Chr == chr, between(Start, left, right))) |>
    select(peak_id = Geneid, Chr, Start, End, log2fc = log2FoldChange, padj, region, locus) |>
    arrange(Chr, Start, End) |>
    filter(region == 31) |> print(n = Inf)

da_unst24_dn24 |>
    select(Geneid, Chr, Start, End, log2FoldChange, padj) |>
    filter(Geneid %in% da_dn24_il24$Geneid) |>
    filter(Chr == "chr7", Start > 50000000, End < 50500000) |>
    arrange(Start, End)


