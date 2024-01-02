# Data handling
library(tidyverse)
library(glue)

# Plotting
library(extrafont)
library(cowplot)
library(patchwork)
library(ggbeeswarm)
library(ggh4x)

#ase_res <- read_tsv("ase_data.tsv", col_types = "ffccciidd")
ase_res <- read_tsv("ase_data_20reads.tsv", col_types = "ffccciidd")

# Plot parameters
stim_colors <- 
    c("Day 0" = "#898e9f",
      "BCR" = "#003967",
      "TLR7" = "#637b31",
      "DN2" = "#a82203")

text_size_big <- 9
text_size_small <- 8

# Plots

# Ref alleles ratio
sample_order <- ase_res |> 
    distinct(sample_id, stim) |> 
    count(sample_id, sort = TRUE) |>
    pull(sample_id)

ref_ratios <- ase_res |>
    mutate(ref_r = ref_count / (ref_count + alt_count),
	   sample_id = factor(sample_id, levels = sample_order)) |>
    select(sample_id, stim, var_id, ref_r)

ref_r_plot <- 
    ggplot(ref_ratios, aes(ref_r)) +
    geom_histogram(aes(fill = stim)) +
    scale_x_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
    scale_fill_manual(values = stim_colors) +
    facet_grid2(sample_id~stim, scales = "free_y", independent = "y") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 11, family = "Arial"),
	  axis.text.y = element_blank(),
	  axis.title.x = element_text(size = 11, family = "Arial"),
	  panel.grid = element_blank(),
	  axis.ticks.x = element_line(linewidth = .5),
	  strip.text.x = element_text(size = 11, family = "Arial"),
	  strip.text.y = element_text(size = 11, angle = 0, family = "Arial"),
	  strip.background = element_rect(fill = "white", color = "white"),
	  legend.position = "none",
	  	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = "Reference allele ratio", y = NULL, fill = "Stim:")

#ggsave("./plots/ref_r.png", ref_r_plot, height = 8, width = 5)
ggsave("./plots/ref_r_20.png", ref_r_plot, height = 8, width = 5)


# Summary stats
summary_1 <- 
    ase_res |>
    group_by(stim) |>
    summarise(n = n_distinct(var_id)) |>
    ungroup()

summary_1_plot <- 
    summary_1 |>
    mutate(labtext = scales::comma(n)) |>
    ggplot(aes(x = stim, y = n)) +
    geom_col(aes(fill = stim), color = "black", linewidth = .5) +
    geom_text(aes(label = labtext),
	      size = 2.5, fontface = "bold", family = "Arial", 
	      vjust = -0.2) +
    scale_fill_manual(values = stim_colors) + 
    theme_minimal() +
    theme(
	  panel.grid.minor.y = element_blank(),
	  panel.grid.major.y = element_blank(),
	  panel.grid.major.x = element_blank(),
	  plot.title = element_text(family = "Arial", size = text_size_big),
	  axis.title.y = element_text(family = "Arial", size = text_size_big),
	  axis.text.x = element_text(family = "Arial", size = text_size_small),
	  axis.text.y = element_blank(),
	  plot.title.position = "plot",
	  plot.background = element_rect(fill = "white", color = "white")
	  ) +
    labs(x = NULL, y = "Number of variants",
	 title = "Number of variants tested per condition\nacross all samples") +
    guides(fill = "none")


summary_2 <- ase_res |>
    group_by(stim, var_id) |>
    slice(which.min(q_value)) |>
    group_by(stim) |>
    summarise("1%" = sum(q_value <= 0.01),
	      "5%" = sum(q_value <= 0.05),
	      "10%" = sum(q_value <= 0.1),
	      "20%" = sum(q_value <= 0.2)) |>
    ungroup() |>
    pivot_longer(-stim, names_to = "FDR", values_to = "n") |>
    mutate(FDR = fct_inorder(FDR))
    
summary_2_plot <- 
    ggplot(summary_2, aes(x = FDR, y = n)) +
    geom_col(aes(fill = stim), position = "dodge", color = "black", linewidth = .5) +
    scale_fill_manual(values = stim_colors) + 
    scale_y_continuous(limits = c(0, max(summary_2$n))) +
    theme_minimal() +
    theme(
	  axis.text = element_text(family = "Arial", size = text_size_small),
	  axis.title = element_text(family = "Arial", size = text_size_big),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white"),
	  plot.title = element_text(family = "Arial", size = text_size_big),
	  plot.title.position = "plot"
	  ) +
    labs(title = "Number of significant variants in any donor\nper FDR threshold",
	 y = "Number of significant variants") +
    guides(fill = "none")

# Genes
summary_3 <- 
    ase_res |>
    filter(!is.na(gene_id)) |>
    select(sample_id, stim, var_id, gene_id, gene_name, q_value) |>
    separate_rows(gene_id, gene_name, sep = "/") |>
    group_by(stim, gene_id) |>
    slice(which.min(q_value)) |> 
    group_by(stim) |>
    summarise("1%" = sum(q_value <= 0.01),
	      "5%" = sum(q_value <= 0.05),
	      "10%" = sum(q_value <= 0.1),
	      "20%" = sum(q_value <= 0.2)) |>
    ungroup() |>
    pivot_longer(-stim, names_to = "FDR", values_to = "n") |>
    mutate(FDR = fct_inorder(FDR))

summary_3_plot <- 
    ggplot(summary_3, aes(x = FDR, y = n)) +
    geom_col(aes(fill = stim), position = "dodge", color = "black", linewidth = .5) +
    scale_fill_manual(values = stim_colors) + 
    scale_y_continuous(limits = c(0, max(summary_3$n))) +
    theme_minimal() +
    theme(
	  axis.text = element_text(family = "Arial", size = text_size_small),
	  axis.title = element_text(family = "Arial", size = text_size_big),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white"),
	  plot.title = element_text(family = "Arial", size = text_size_big),
	  plot.title.position = "plot"
	  ) +
    labs(title = "Number of unique genes with significant ASE variants\nper FDR threshold",
	 y = "Number of genes") +
    guides(fill = "none")


# Summary per person
summary_4 <- 
    ase_res |>
    count(stim, sample_id)

summary_4_plot <- 
    ggplot(summary_4, aes(x = stim, y = n)) +
    geom_quasirandom(aes(color = stim), method = "smiley", width = .25, alpha = .75) +
    geom_boxplot(aes(fill = stim), linewidth = .25, alpha = .5, width = .5, outlier.color = NA) +
    scale_y_continuous(limits = c(0, max(summary_4$n)), label = scales::comma) +
    scale_color_manual(values = stim_colors) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(
	  axis.text = element_text(family = "Arial", size = text_size_small),
	  axis.title = element_text(family = "Arial", size = text_size_big),
	  panel.grid.major.x = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white"),
	  plot.title = element_text(family = "Arial", size = text_size_big),
	  plot.title.position = "plot"
	  ) +
    labs(x = NULL, y = "Number of variants",
	 title = "Total number of variants tested\nper sample") +
    guides(fill = "none", color = "none")

summary_5 <- 
    ase_res |>
    group_by(stim, sample_id) |>
    summarise("1%" = sum(q_value <= 0.01),
	      "5%" = sum(q_value <= 0.05),
	      "10%" = sum(q_value <= 0.1),
	      "20%" = sum(q_value <= 0.2)) |>
    ungroup() |>
    pivot_longer(-c(stim, sample_id), names_to = "FDR", values_to = "n") |>
    mutate(FDR = fct_inorder(FDR))

summary_5_plot <- 
    ggplot(summary_5, aes(x = stim, y = n)) +
    geom_quasirandom(aes(color = stim), 
		     method = "smiley", width = .2, size = .75, alpha = .75) +
    geom_boxplot(aes(fill = stim), linewidth = .25, alpha = .5, width = .5, outlier.color = NA) +
    scale_y_continuous(limits = c(0, max(summary_5$n)), 
		       breaks = seq(0, 1400, 200),
		       label = scales::comma) +
    scale_color_manual(values = stim_colors) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~FDR, nrow = 1, strip.position = "bottom") +
    theme_minimal() +
    theme(
	  axis.text.x = element_blank(),
	  axis.text.y = element_text(family = "Arial", size = text_size_small),
	  axis.title.x = element_text(family = "Arial", size = text_size_big, margin = margin(t = -10)),
	  axis.title.y = element_text(family = "Arial", size = text_size_big),
	  strip.text = element_text(family = "Arial", size = text_size_small),
	  plot.background = element_rect(color = "white", fill = "white"),
	  panel.spacing = unit(0.2, "lines"),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.title = element_text(family = "Arial", size = text_size_big),
	  plot.title.position = "plot"
	  ) +
    labs(x = "FDR", y = "Number of variants",
	 title = "Number of significant variants per sample\nper FDR threshold") +
    guides(color = "none", fill = "none")


summary_6 <- ase_res |>
    filter(q_value <= 0.2) |>
    separate_rows(gene_id, gene_name, sep = "/") |>
    group_by(stim, sample_id, gene_id) |>
    summarise("0.01" = any(q_value <= 0.01),
	      "0.05" = any(q_value <= 0.05),
	      "0.1" = any(q_value <= 0.1),
	      "0.2" = any(q_value <= .2)) |>
    ungroup() |>
    pivot_longer(-c(stim, sample_id, gene_id), names_to = "FDR") |>
    group_by(stim, sample_id, FDR) |>
    summarise(n = sum(value)) |>
    ungroup()

summary_6_plot <- 
    ggplot(summary_6, aes(x = stim, y = n)) +
    geom_quasirandom(aes(color = stim), 
		     method = "smiley", width = .2, size = .75, alpha = .75) +
    geom_boxplot(aes(fill = stim), linewidth = .25, alpha = .5, width = .5, outlier.color = NA) +
    scale_y_continuous(limits = c(0, max(summary_6$n)), 
		       breaks = seq(0, 800, 200)) +
    scale_color_manual(values = stim_colors) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~FDR, nrow = 1, strip.position = "bottom") +
    theme_minimal() +
    theme(
	  axis.text.x = element_blank(),
	  axis.text.y = element_text(family = "Arial", size = text_size_small),
	  axis.title.x = element_text(family = "Arial", size = text_size_big, margin = margin(t = -10)),
	  axis.title.y = element_text(family = "Arial", size = text_size_big),
	  strip.text = element_text(family = "Arial", size = text_size_small),
	  plot.background = element_rect(color = "white", fill = "white"),
	  panel.spacing = unit(0.2, "lines"),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.title = element_text(family = "Arial", size = text_size_big),
	  plot.title.position = "plot"
	  ) +
    labs(x = "FDR", y = "Number of variants",
	 title = "Number of genes with significant ASE variants\nper sample per FDR threshold") +
    guides(color = "none", fill = "none")


plot_out <- 
    summary_1_plot + summary_2_plot + summary_3_plot +
    summary_4_plot + summary_5_plot + summary_6_plot + 
    plot_layout(ncol = 3, widths = c(1, 1.5, 1.5))

#ggsave("./plots/summary.png", plot_out, height = 5.5, width = 8.5, dpi = 600)
ggsave("./plots/summary_20.png", plot_out, height = 5.5, width = 8.5, dpi = 600)


## QQ plot
make_qq_df <- function(ps) {

    n <- length(ps)

    tibble(observed = -log10(sort(ps)),
	   expected = -log10(ppoints(n)))
}

qq_df <- ase_res |>
    unite("id", c(sample_id, stim), sep = "_") |>
    mutate(p_value = ifelse(p_value < 5e-324, 5e-324, p_value)) |>
    {function(x) split(x, x$id)}() |>
    map(~pull(., p_value)) |>
    map_df(make_qq_df, .id = "id") |>
    separate(id, c("sample_id", "stim"), sep = "_") |>
    mutate(sample_id = factor(sample_id, levels = levels(ase_res$sample_id)),
	   stim = factor(stim, levels = levels(ase_res$stim)))

qq_plot <-
    ggplot(qq_df) +
    geom_point(aes(x = expected, y = observed, color = stim), size = 1, alpha = .5) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = stim_colors) +
    facet_grid(sample_id~stim) +
    theme_minimal() +
    theme(
	  axis.text = element_text(family = "Arial", size = 11),
	  axis.title = element_text(family = "Arial", size = 14),
	  strip.text = element_text(family = "Arial", size = 14),
	  panel.grid = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = expression(paste("Expected -log"[10], plain(P))),
	 y = expression(paste("Observed -log"[10], plain(P)))) +
    guides(color = "none")

#ggsave("./plots/qq.png", qq_plot, height = 12, width = 5, dpi = 600)
ggsave("./plots/qq_20.png", qq_plot, height = 12, width = 5, dpi = 600)

p_hist <- 
    ggplot(ase_res, aes(x = p_value)) +
    geom_histogram(aes(fill = stim), bins = 20, linewidth = 0) +
    scale_fill_manual(values = stim_colors) +
    scale_x_continuous(breaks = c(0, 1)) +
    facet_grid2(sample_id~stim, scales = "free_y", independent = "y") +
    theme_minimal() +
    theme(
	  axis.text.x = element_text(family = "Arial", size = 11),
	  axis.text.y = element_blank(),
	  axis.title = element_text(family = "Arial", size = 14),
	  strip.text = element_text(family = "Arial", size = 14),
	  panel.grid = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  plot.background = element_rect(color = "white", fill = "white")) +
    guides(fill = "none") +
    labs(x = "p-value", y = NULL)

#ggsave("./plots/pvalue_hist.png", p_hist, height = 12, width = 5, dpi = 600)
ggsave("./plots/pvalue_hist_20.png", p_hist, height = 12, width = 5, dpi = 600)

# P-value ~ coverage
testplot <- 
    ase_res |>
    filter(sample_id == first(sample_id)) |>
    mutate(total = ref_count + alt_count) |>
    select(sample_id, stim, p_value, total) |>
    mutate(bin = cut_interval(p_value, n = 20)) |>
    arrange(bin) |>
    ggplot(aes(x = bin, y = log2(total))) +
	geom_violin(aes(fill = stim)) +
	facet_wrap(~stim, ncol = 2) +
	#facet_grid2(sample_id~stim, scales = "free_y", independent = "y") +
	scale_fill_manual(values = stim_colors) +
	theme_bw() +
	theme(
	      axis.text.x = element_text(family = "Arial", size = 11, angle = 45, hjust = 1, vjust = 1),
	      axis.text.y = element_text(family = "Arial", size = 11),
	      axis.title = element_text(family = "Arial", size = 14),
	      strip.text = element_text(family = "Arial", size = 14),
	      panel.grid = element_blank(),
	      panel.border = element_blank(),
	      strip.text.y = element_text(angle = 0),
	      plot.background = element_rect(color = "white", fill = "white")) +
	guides(fill = "none") +
	labs(x = "P-value bins", y = "Log2 (number of reads)")

#ggsave("./plots/testplot.png", testplot, width = 8, height = 8)
ggsave("./plots/testplot_20.png", testplot, width = 8, height = 8)




# Number of unique reads
temp_work <- system("echo $TEMP_WORK", intern = TRUE)

log_df <- 
    ase_res |>
    distinct(sample_id, stim) |>
    mutate(stim2 = recode(stim, "Day 0" = "unstday0")) |>
    mutate(logfile = file.path(temp_work, glue("/bam/fixase/{sample_id}_{stim2}_Log.final.out")),
	   uniq = map_int(logfile, ~read_tsv(., col_names = FALSE) |>
			  filter(grepl("^Uniquely mapped reads number", X1)) |>
			  pull(X2) |>
			  as.integer())) |>
    select(sample_id, stim, uniq)

# Correlation between number of significant hits and total number of reads
ase_by_coverage_plot <- 
    left_join(log_df, filter(summary_5, FDR == "10%"), join_by(sample_id, stim)) |>
    select(sample_id, stim, n_reads = uniq, n_ase = n) |>
    ggplot(aes(x = n_ase, y = n_reads)) +
	geom_point(aes(color = stim), size = 4, alpha = .7) +
	scale_color_manual(values = stim_colors) +
	scale_x_continuous(limits = c(0, 1000)) +
	scale_y_continuous(labels = function(x) x/1e6L) +
	theme_minimal() +
	theme(
	      axis.text = element_text(family = "Arial", size = 11),
	      axis.title = element_text(family = "Arial", size = 14),
	      panel.grid.major = element_line(color = "grey96"),
	      panel.grid.minor = element_blank(),
	      plot.background = element_rect(color = "white", fill = "white"),
	      plot.title = element_text(family = "Arial", size = 14),
	      plot.title.position = "plot") +
	labs(x = "Number of ASE variants", 
	     y = "Uniquely mapped reads (Millions)",
	     color = "Stim:",
	     title = "Number of ASE variants at 10% FDR\nvs. number of uniquely mapped reads")

#ggsave("./plots/ase_by_coverage.png", ase_by_coverage_plot, width = 6, height = 5)
ggsave("./plots/ase_by_coverage_20.png", ase_by_coverage_plot, width = 6, height = 5)

# Number of SLE genes with ASE per donor
# When there are more than 1 gene per locus, I followed the rule below:
# - For Langefeld et al, I think it makes more sense to keep the first gene most of the times;
# - For Bentham et al, I think it makes sense to keep both.

langefeld <- 
    "./imbalance_screening/data/langefeld_hits.tsv" |>
    read_tsv() |>
    distinct(gene) |>
    pull(gene) |>
    str_split("-") |>
    map_chr(function(x) ifelse(any(grepl("IKZF", x)), x[grepl("IKZ", x)], x[1]))

bentham <- "../../../sle_variants/paper_data/bentham_tab1.tsv" |>
    read_tsv() |>
    select(locus) |>
    filter(!grepl("^MHC", locus)) |>
    separate_rows(locus, sep = ", ") |>
    mutate(locus = recode(locus, "CXorf21" = "TASL")) |>
    pull(locus)

gwas_genes <- 
    c(langefeld, bentham) |>
    unique() |>
    c("STAT1") |>
    sort()
    
ase_load <- ase_res |>
    filter(q_value <= 0.05) |>
    separate_rows(gene_id, gene_name, sep = "/") |>
    filter(gene_name %in% gwas_genes) |>
    count(sample_id, stim, gene_name) |>
    mutate(gene_name = factor(gene_name, levels = gwas_genes)) |>
    complete(sample_id, stim, gene_name, fill = list(n = 0)) |>
    arrange(sample_id, stim, gene_name) |>
    group_by(gene_name) |>
    filter(any(n > 0)) |>
    ungroup()


plot_tiles <- function(stim_name) {

    plot_data <- filter(ase_load, stim == stim_name)

    ggplot(plot_data, aes(x = sample_id, y = gene_name)) +
    geom_tile(aes(fill = n), alpha = .5) +
    geom_text(data = plot_data |> filter(n > 0),
	      aes(label = n), family = "Arial") +
    scale_fill_gradient(low = "white", high = stim_colors[[stim_name]]) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, family = "Arial", angle = 45, hjust = 1, vjust = 1),
	  axis.text.y = element_text(size = 12, family = "Arial"),
	  panel.border = element_blank(),
	  plot.title = element_text(size = 16, family = "Arial", hjust = .5),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = NULL, y = NULL, title = stim_name) +
    guides(fill = "none")
}

tiles <- 
    plot_tiles("Day 0") + 
    plot_tiles("TLR7") + 
    plot_tiles("BCR") + 
    plot_tiles("DN2") + 
    plot_layout(nrow = 1)

ggsave("./plots/tiles.png", tiles, width = 24, height = 7, dpi = 600)

################################################################################
# Replication
compute_replication <- function(sample_1, sample_2) {
    
    temp1 <- 
	ref_ratios |>
	filter(sample_id == sample_1, q_value <= 0.05) |>
	select(sample_id:ref_r)

    ref_ratios |>
	filter(sample_id == sample_2) |>
	inner_join(temp1, join_by(stim, var_id)) |>
	group_by(stim) |>
	summarise(r = cor(ref_r.x, ref_r.y)) |>
	ungroup()
}

ref_ratios <- 
    ase_res |> 
    mutate(ref_r = refCount/total) |>
    select(sample_id, stim, var_id, ref_r, q_value)

cor_df <- 
    ref_ratios |>
    distinct(sample_id) |>
    expand_grid(sample_1 = sample_id, sample_2 = sample_id) |>
    distinct(sample_1, sample_2) |>
    mutate(data = map2(sample_1, sample_2, compute_replication)) |>
    unnest(cols = c(data))


plot_corr <- function(stim_i) {

    lim <- round(range(cor_df$r), 1)

    ggplot(cor_df |> filter(stim == stim_i), 
       aes(x = sample_1, y = sample_2)) +
    geom_tile(aes(fill = r), alpha = .8) +
    scale_fill_gradient(low = "white", high = stim_colors[stim_i]) +
    geom_text(aes(label = round(r, 2)),
	      size = 4, fontface = "bold", family = "Arial") +
    scale_x_discrete(labels = function(x) sub("^(\\d+)\\.(\\d)$", "\\1 rep #\\2", x))  +
    scale_y_discrete(labels = function(x) sub("^(\\d+)\\.(\\d)$", "\\1 rep #\\2", x))  +
    theme_minimal() +
    theme(
	  axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, family = "Arial"),
	  axis.text.y = element_text(size = 12, family = "Arial"),
	  legend.text = element_text(size = 12, family = "Arial"),  
	  legend.title = element_text(size = 12, family = "Arial"),  
	  plot.title = element_text(size = 16, hjust = .5, family = "Arial", face = "bold"),
	  panel.grid = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = NULL, y = NULL, title = stim_i) +
    guides(fill = "none")
}


corr_plot <- 
    plot_corr("Day 0") + plot_corr("TLR7") + plot_corr("BCR") + plot_corr("DN2") +
    plot_layout(ncol = 2) +
    plot_annotation(title = "Spearman correlation of reference allele ratios among samples at variants with significant ASE at 5% FDR.",
		    theme = theme(plot.title = element_text(size = 24, family = "Arial")))

ggsave("./plots/reps_corr.png", corr_plot, height = 18, width = 18, dpi = 300)

