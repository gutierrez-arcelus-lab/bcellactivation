# Data handling
library(tidyverse)
library(tidytext)
library(glue)

# Plotting
library(extrafont)
library(cowplot)
library(patchwork)
library(ggbeeswarm)
library(ggh4x)


if (!file.exists("plots")) dir.create("plots")

# Plot parameters
stim_colors <- 
    c("Day 0" = "#898e9f",
      "BCR" = "#003967",
      "TLR7" = "#637b31",
      "DN2" = "#a82203")

ase_res <- 
    read_tsv("ase_data.tsv", col_types = "ffcddfccdd") |>
    mutate(fdr = p.adjust(p_value, method = "fdr")) |>
    select(-q_value)

# Plots
# Ref alleles ratio
ref_ratios <- ase_res |>
    mutate(ref_r = ref_count / (ref_count + alt_count)) |>
    select(sample_id, stim, variant_id, ref_r)

ref_ratio_plot <- 
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

ggsave("./plots/ref_ratios.png", ref_ratio_plot, height = 8, width = 5)

# Summary stats
summary_1 <- 
    ase_res |>
    group_by(stim) |>
    summarise(n = n_distinct(variant_id)) |>
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
	  plot.title = element_text(family = "Arial", size = 11),
	  axis.title.y = element_text(family = "Arial", size = 11),
	  axis.text.x = element_text(family = "Arial", size = 9),
	  axis.text.y = element_blank(),
	  plot.title.position = "plot",
	  plot.background = element_rect(fill = "white", color = "white")
	  ) +
    labs(x = NULL, y = NULL,
	 title = "Number of variants tested per\ncondition across all samples") +
    guides(fill = "none")


summary_2 <- ase_res |>
    group_by(stim, variant_id) |>
    slice(which.min(fdr)) |>
    group_by(stim) |>
    summarise("1%" = sum(fdr <= 0.01),
	      "5%" = sum(fdr <= 0.05),
	      "10%" = sum(fdr <= 0.1),
	      "20%" = sum(fdr <= 0.2)) |>
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
	  axis.text = element_text(family = "Arial", size = 9),
	  axis.title = element_text(family = "Arial", size = 11),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white"),
	  plot.title = element_text(family = "Arial", size = 11),
	  plot.title.position = "plot"
	  ) +
    labs(title = "Number of significant variants in any donor",
	 y = NULL) +
    guides(fill = "none")

# Genes
summary_3 <- 
    ase_res |>
    filter(annot == "exon") |>
    select(sample_id, stim, variant_id, gene_id, gene_name, fdr) |>
    separate_rows(gene_id, gene_name, sep = "/") |>
    group_by(stim, gene_id) |>
    slice(which.min(fdr)) |> 
    group_by(stim) |>
    summarise("1%" = sum(fdr <= 0.01),
	      "5%" = sum(fdr <= 0.05),
	      "10%" = sum(fdr <= 0.1),
	      "20%" = sum(fdr <= 0.2)) |>
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
	  axis.text = element_text(family = "Arial", size = 9),
	  axis.title = element_text(family = "Arial", size = 11),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white"),
	  plot.title = element_text(family = "Arial", size = 11),
	  plot.title.position = "plot"
	  ) +
    labs(title = "Number of unique genes with\nsignificant ASE variants",
	 y = NULL) +
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
	  axis.text = element_text(family = "Arial", size = 9),
	  axis.title = element_text(family = "Arial", size = 11),
	  panel.grid.major.x = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white"),
	  plot.title = element_text(family = "Arial", size = 11),
	  plot.title.position = "plot"
	  ) +
    labs(x = NULL, y = NULL,
	 title = "Total number of variants tested\nper sample") +
    guides(fill = "none", color = "none")

summary_5 <- 
    ase_res |>
    group_by(stim, sample_id) |>
    summarise("1%" = sum(fdr <= 0.01),
	      "5%" = sum(fdr <= 0.05),
	      "10%" = sum(fdr <= 0.1),
	      "20%" = sum(fdr <= 0.2)) |>
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
	  axis.text.y = element_text(family = "Arial", size = 9),
	  axis.title.x = element_text(family = "Arial", size = 11, margin = margin(t = -10)),
	  axis.title.y = element_text(family = "Arial", size = 11),
	  strip.text = element_text(family = "Arial", size = 9),
	  plot.background = element_rect(color = "white", fill = "white"),
	  panel.spacing = unit(0.2, "lines"),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.title = element_text(family = "Arial", size = 11),
	  plot.title.position = "plot"
	  ) +
    labs(x = "FDR", y = NULL,
	 title = "Number of significant variants per sample") +
    guides(color = "none", fill = "none")


summary_6 <- ase_res |>
    filter(fdr <= 0.2) |>
    separate_rows(gene_id, gene_name, sep = "/") |>
    group_by(stim, sample_id, gene_id) |>
    summarise("0.01" = any(fdr <= 0.01),
	      "0.05" = any(fdr <= 0.05),
	      "0.1" = any(fdr <= 0.1),
	      "0.2" = any(fdr <= .2)) |>
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
	  axis.text.y = element_text(family = "Arial", size = 9),
	  axis.title.x = element_text(family = "Arial", size = 11, margin = margin(t = -10)),
	  axis.title.y = element_text(family = "Arial", size = 11),
	  strip.text = element_text(family = "Arial", size = 9),
	  plot.background = element_rect(color = "white", fill = "white"),
	  panel.spacing = unit(0.2, "lines"),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.title = element_text(family = "Arial", size = 11),
	  plot.title.position = "plot"
	  ) +
    labs(x = "FDR", y = NULL,
	 title = "Number of genes with significant\nASE variants per sample") +
    guides(color = "none", fill = "none")


plot_out <- 
    summary_1_plot + summary_2_plot + summary_3_plot +
    summary_4_plot + summary_5_plot + summary_6_plot + 
    plot_layout(ncol = 3, widths = c(1, 1.5, 1.5))

ggsave("./plots/summary.png", plot_out, height = 5.5, width = 8.5, dpi = 600)


## QQ plot
make_qq_df <- function(ps) {

    n <- length(ps)

    tibble(observed = -log10(sort(ps)),
	   expected = -log10(ppoints(n)))
}

qq_df <- 
    ase_res |>
    mutate(p_value = ifelse(p_value < 5e-324, 5e-324, p_value)) |>
    {function(x) split(x, list(x$sample_id, x$stim))}() |>
    map(~pull(., p_value)) |>
    map_df(make_qq_df, .id = "id") |>
    separate(id, c("sample_id", "stim"), sep = "\\.") |>
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

ggsave("./plots/qq.png", qq_plot, height = 12, width = 5, dpi = 600)

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

ggsave("./plots/pvalue_hist.png", p_hist, height = 12, width = 5, dpi = 600)


# Correlation between number of significant hits and total number of reads
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
	mutate(reads = tolower(reads),
	       reads = factor(reads, levels = c("unique", "multi-mapping", "unmapped", "chimeric"))) |>
	arrange(reads)
}

temp_work <- system("echo $TEMP_WORK", intern = TRUE)

star_log_df <- 
    ase_res |>
    distinct(sample_id, stim) |>
    mutate(stim2 = recode(stim, "Day 0" = "unstday0")) |>
    mutate(logfile = file.path(temp_work, glue("/bam/ase_v2/{sample_id}_{stim2}/Log.final.out")),
	   data = map(logfile, read_star_log)) |>
    select(sample_id, stim, data) |>
    unnest(cols = data)

ase_by_coverage_plot <- 
    left_join(filter(star_log_df, reads == "unique"),
	      filter(summary_5, FDR == "10%"), 
	      join_by(sample_id, stim)) |>
    ggplot(aes(x = n, y = number_of_reads)) +
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

ggsave("./plots/ase_by_coverage.png", ase_by_coverage_plot, width = 6, height = 5)

###############################################################################
# Plot glm

## GLM results
res_norand_df <- 
    read_tsv("./results_glm/glm_res_df_annotated.tsv", col_types = "cfccddccc") |>
    filter(grepl("^Day 0-", stim)) |>
    mutate(stim = str_remove(stim, "Day 0-"),
	   stim = factor(stim, levels = names(stim_colors)[-1]),
	   fdr = p.adjust(p, method = "fdr"),
	   gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    select(donor_id, replic, stim, variant_id, annot, gene_id, gene_name, beta_stim, p, fdr)

## Summary
summary_1_dyn <- 
    res_norand_df |>
    group_by(stim) |>
    summarise(n = n_distinct(variant_id)) |>
    ungroup()

summary_1_dyn_plot <- 
    summary_1_dyn |>
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
	  plot.title = element_text(family = "Arial", size = 11),
	  axis.title.y = element_text(family = "Arial", size = 11),
	  axis.text.x = element_text(family = "Arial", size = 9),
	  axis.text.y = element_blank(),
	  plot.title.position = "plot",
	  plot.background = element_rect(fill = "white", color = "white")
	  ) +
    labs(x = NULL, y = NULL,
	 title = "Number of variants tested per\ncondition across all samples") +
    guides(fill = "none")


summary_2_dyn <- 
    res_norand_df |>
    group_by(stim, variant_id) |>
    slice(which.min(fdr)) |>
    group_by(stim) |>
    summarise("1%" = sum(fdr <= 0.01),
	      "5%" = sum(fdr <= 0.05),
	      "10%" = sum(fdr <= 0.1),
	      "20%" = sum(fdr <= 0.2)) |>
    ungroup() |>
    pivot_longer(-stim, names_to = "FDR", values_to = "n") |>
    mutate(FDR = fct_inorder(FDR))
    
summary_2_dyn_plot <- 
    ggplot(summary_2_dyn, aes(x = FDR, y = n)) +
    geom_col(aes(fill = stim), position = "dodge", color = "black", linewidth = .5) +
    scale_fill_manual(values = stim_colors) + 
    scale_y_continuous(limits = c(0, max(summary_2_dyn$n))) +
    theme_minimal() +
    theme(
	  axis.text = element_text(family = "Arial", size = 9),
	  axis.title = element_text(family = "Arial", size = 11),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white"),
	  plot.title = element_text(family = "Arial", size = 11),
	  plot.title.position = "plot"
	  ) +
    labs(title = "Number of significant variants in any donor",
	 y = NULL) +
    guides(fill = "none")


# Genes
summary_3_dyn <- 
    res_norand_df |>
    filter(annot == "exon") |>
    unite("sample_id", c(donor_id, replic), sep = "_") |>
    select(sample_id, stim, variant_id, gene_id, gene_name, fdr) |>
    separate_rows(gene_id, gene_name, sep = "/") |>
    group_by(stim, gene_id) |>
    slice(which.min(fdr)) |> 
    group_by(stim) |>
    summarise("1%" = sum(fdr <= 0.01),
	      "5%" = sum(fdr <= 0.05),
	      "10%" = sum(fdr <= 0.1),
	      "20%" = sum(fdr <= 0.2)) |>
    ungroup() |>
    pivot_longer(-stim, names_to = "FDR", values_to = "n") |>
    mutate(FDR = fct_inorder(FDR))

summary_3_dyn_plot <- 
    ggplot(summary_3_dyn, aes(x = FDR, y = n)) +
    geom_col(aes(fill = stim), position = "dodge", color = "black", linewidth = .5) +
    scale_fill_manual(values = stim_colors) + 
    scale_y_continuous(limits = c(0, max(summary_3_dyn$n))) +
    theme_minimal() +
    theme(
	  axis.text = element_text(family = "Arial", size = 9),
	  axis.title = element_text(family = "Arial", size = 11),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white"),
	  plot.title = element_text(family = "Arial", size = 11),
	  plot.title.position = "plot"
	  ) +
    labs(title = "Number of unique genes with\nsignificant ASE variants",
	 y = NULL) +
    guides(fill = "none")

# Summary per person
summary_4_dyn <- 
    res_norand_df |>
    unite("sample_id", c(donor_id, replic), sep = "_") |>
    count(stim, sample_id)

summary_4_dyn_plot <- 
    ggplot(summary_4_dyn, aes(x = stim, y = n)) +
    geom_quasirandom(aes(color = stim), method = "smiley", width = .25, alpha = .75) +
    geom_boxplot(aes(fill = stim), linewidth = .25, alpha = .5, width = .5, outlier.color = NA) +
    scale_y_continuous(limits = c(0, max(summary_4_dyn$n)), label = scales::comma) +
    scale_color_manual(values = stim_colors) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(
	  axis.text = element_text(family = "Arial", size = 9),
	  axis.title = element_text(family = "Arial", size = 11),
	  panel.grid.major.x = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white"),
	  plot.title = element_text(family = "Arial", size = 11),
	  plot.title.position = "plot"
	  ) +
    labs(x = NULL, y = NULL,
	 title = "Total number of variants tested\nper sample") +
    guides(fill = "none", color = "none")

summary_5_dyn <- 
    res_norand_df |>
    unite("sample_id", c(donor_id, replic), sep = "_") |>
    group_by(stim, sample_id) |>
    summarise("1%" = sum(fdr <= 0.01),
	      "5%" = sum(fdr <= 0.05),
	      "10%" = sum(fdr <= 0.1),
	      "20%" = sum(fdr <= 0.2)) |>
    ungroup() |>
    pivot_longer(-c(stim, sample_id), names_to = "FDR", values_to = "n") |>
    mutate(FDR = fct_inorder(FDR))

summary_5_dyn_plot <- 
    ggplot(summary_5_dyn, aes(x = stim, y = n)) +
    geom_quasirandom(aes(color = stim), 
		     method = "smiley", width = .2, size = .75, alpha = .75) +
    geom_boxplot(aes(fill = stim), linewidth = .25, alpha = .5, width = .5, outlier.color = NA) +
    scale_y_continuous(limits = c(0, max(summary_5_dyn$n)), 
		       label = scales::comma) +
    scale_color_manual(values = stim_colors) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~FDR, nrow = 1, strip.position = "bottom") +
    theme_minimal() +
    theme(
	  axis.text.x = element_blank(),
	  axis.text.y = element_text(family = "Arial", size = 9),
	  axis.title.x = element_text(family = "Arial", size = 11, margin = margin(t = -10)),
	  axis.title.y = element_text(family = "Arial", size = 11),
	  strip.text = element_text(family = "Arial", size = 9),
	  plot.background = element_rect(color = "white", fill = "white"),
	  panel.spacing = unit(0.2, "lines"),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.title = element_text(family = "Arial", size = 11),
	  plot.title.position = "plot"
	  ) +
    labs(x = "FDR", y = NULL,
	 title = "Number of significant variants per sample") +
    guides(color = "none", fill = "none")


summary_6_dyn <- 
    res_norand_df |>
    unite("sample_id", c(donor_id, replic), sep = "_") |>
    filter(fdr <= 0.2) |>
    separate_rows(gene_id, gene_name, sep = "/") |>
    group_by(stim, sample_id, gene_id) |>
    summarise("0.01" = any(fdr <= 0.01),
	      "0.05" = any(fdr <= 0.05),
	      "0.1" = any(fdr <= 0.1),
	      "0.2" = any(fdr <= .2)) |>
    ungroup() |>
    pivot_longer(-c(stim, sample_id, gene_id), names_to = "FDR") |>
    group_by(stim, sample_id, FDR) |>
    summarise(n = sum(value)) |>
    ungroup()

summary_6_dyn_plot <- 
    ggplot(summary_6_dyn, aes(x = stim, y = n)) +
    geom_quasirandom(aes(color = stim), 
		     method = "smiley", width = .2, size = .75, alpha = .75) +
    geom_boxplot(aes(fill = stim), linewidth = .25, alpha = .5, width = .5, outlier.color = NA) +
    scale_y_continuous(limits = c(0, max(summary_6_dyn$n))) +
    scale_color_manual(values = stim_colors) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~FDR, nrow = 1, strip.position = "bottom") +
    theme_minimal() +
    theme(
	  axis.text.x = element_blank(),
	  axis.text.y = element_text(family = "Arial", size = 9),
	  axis.title.x = element_text(family = "Arial", size = 11, margin = margin(t = -10)),
	  axis.title.y = element_text(family = "Arial", size = 11),
	  strip.text = element_text(family = "Arial", size = 9),
	  plot.background = element_rect(color = "white", fill = "white"),
	  panel.spacing = unit(0.2, "lines"),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.title = element_text(family = "Arial", size = 11),
	  plot.title.position = "plot"
	  ) +
    labs(x = "FDR", y = NULL,
	 title = "Number of genes with significant\nASE variants per sample") +
    guides(color = "none", fill = "none")


plot_out_dyn <- 
    summary_1_dyn_plot + summary_2_dyn_plot + summary_3_dyn_plot +
    summary_4_dyn_plot + summary_5_dyn_plot + summary_6_dyn_plot + 
    plot_layout(ncol = 3, widths = c(1, 1.5, 1.5))

ggsave("./plots/summary_dyn.png", plot_out_dyn, height = 5.5, width = 8.5, dpi = 600)



## QQ-plot
qq_norand_df <- 
    res_norand_df |>
    filter(!is.na(p)) |>
    group_split(donor_id, replic, stim) |>
    map_df(~arrange(., p) |> 
	   mutate(observed = -log10(p), 
		  expected = -log10(ppoints(n())))) |>
    unite("sample_id", c(donor_id, replic), sep = "_")

qq_plot <- 
    ggplot(qq_norand_df) +
    geom_point(aes(x = expected, y = observed, color = stim),
	       size = 1, alpha = .5) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = stim_colors) +
    facet_grid(sample_id~stim, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "none",
	  panel.grid.minor = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = expression(paste("Expected -log"[10], plain(P))),
	 y = expression(paste("Observed -log"[10], plain(P))))

ggsave("./plots/qq_norand_glm.png", qq_plot, height = 12, width = 4)

# Plot betas between technical replicates
get_rep_betas <- function(d, s, r1, r2) {
    
    data_r1 <- 
	res_norand_df |>
	filter(donor_id == d, replic == r1, stim == s,
	       fdr <= 0.1) |>
	select(donor_id, stim, variant_id, beta_stim)

    data_r2 <-
	res_norand_df |>
	filter(donor_id == d, replic == r2, stim == s) |>
	inner_join(data_r1, join_by(donor_id, stim, variant_id)) |>
	select(variant_id, beta_stim.x, beta_stim.y, p)

    data_r2
}

replic_data <- 
    res_norand_df |>
    distinct(donor_id = donor_id, replic, stim) |>
    group_by(donor_id, stim) |>
    filter(n_distinct(replic) > 1) |>
    ungroup()

replic_spec <-     
    left_join(replic_data, replic_data, 
	  join_by(donor_id, stim), 
	  relationship = "many-to-many") |>
    filter(replic.x != replic.y) |>
    select(donor_id, stim, replic.x, replic.y) |>
    arrange(donor_id, stim, replic.x, replic.y) |>
    mutate(data = pmap(list(donor_id, stim, replic.x, replic.y), get_rep_betas)) |>
    unnest(data)

plot_data <- 
    replic_spec |>
    mutate(lab = sprintf("Donor #%s\n Rep %s vs. Rep %s", donor_id, replic.x, replic.y))

cor_data <- 
    plot_data |>
    group_by(lab) |>
    mutate(ymax = max(beta_stim.y)) |>
    group_by(stim) |>
    mutate(xmin = min(beta_stim.x)) |>
    group_by(lab, donor_id, stim, replic.x, replic.y, xmin, ymax) |>
    summarise(r = cor(x = beta_stim.x, y = beta_stim.y, method = "spearman")) |>
    ungroup() |>
    mutate(cor_lab = paste("\u03C1 =", round(r, 2)))

p_beta_reps <- 
    ggplot(plot_data, aes(x = beta_stim.x, y = beta_stim.y)) +
    geom_abline() +
    geom_vline(xintercept = 0, color = "grey90") +
    geom_hline(yintercept = 0, color = "grey90") +
    geom_point(aes(color = stim), alpha = .75) +
    scale_color_manual(values = stim_colors) +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    geom_text(data = cor_data, aes(x = xmin, y = ymax + .25, label = cor_lab),
	      size = 3.5, family = "Arial", hjust = "inward") +
    facet_grid(lab~stim, scales = "free") +
    coord_cartesian(clip = 'off') +
    theme_minimal() +
    theme(legend.position = "none",
	  axis.title = element_blank(),
	  strip.text.x = element_text(size = 12, family = "arial", face = "bold"),
	  strip.text.y = element_text(size = 12, family = "arial", angle = 0),
	  panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/betas_replic.png", p_beta_reps, height = 8.5, width = 5)

# test
min_cov_df <- 
    ase_res |>
    separate(sample_id, c("donor_id", "replic"), sep = "_") |>
    mutate(replic = recode(replic, "1" = "A", "2" = "B", "3" = "C"),
	   total = ref_count + alt_count) |>
    select(donor_id, replic, stim, variant_id, total) |>
    pivot_wider(names_from = stim, values_from = total) |>
    pivot_longer(BCR:DN2, names_to = "stim", values_to = "total") |>
    drop_na(total) |>
    mutate(min_cov = pmin(`Day 0`, total)) |>
    select(donor_id, replic, stim, variant_id, min_cov)

plot_data_cov <- 
    left_join(plot_data, min_cov_df, join_by(donor_id, stim, variant_id, replic.x == replic)) |>
    left_join(min_cov_df, join_by(donor_id, stim, variant_id, replic.y == replic)) |>
    mutate(min_cov = pmin(min_cov.x, min_cov.y)) |>
    select(-min_cov.x, -min_cov.y) |>
    mutate(stim = factor(stim, levels = c("BCR", "TLR7", "DN2")))

p_beta_reps_cov <- 
    ggplot(plot_data_cov, 
	   aes(x = beta_stim.x, 
	       y = beta_stim.y)) +
    geom_abline() +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_point(aes(color = log2(min_cov)), size = .7) +
    scale_color_viridis_c("log2 min counts:", 
			  guide = guide_colorbar(barheight = .5)) +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    geom_text(data = cor_data, aes(x = xmin, y = ymax + .1, label = cor_lab),
	      size = 3, family = "Arial", hjust = "inward") +
    facet_grid(lab~stim, scales = "free") +
    theme_bw() +
    theme(legend.position = 'top',
	  axis.title = element_blank(),
	  strip.text.x = element_text(size = 12, family = "Arial", face = "bold"),
	  strip.text.y = element_text(size = 12, family = "Arial", angle = 0),
	  panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/betas_replic_cov.png", p_beta_reps_cov, height = 8.5, width = 5)



# Plot cases of significant effects of stim
fdr01 <- 
    res_norand_df |>
    group_by(variant_id) |>
    filter(any(fdr <= 0.01)) |>
    ungroup() |>
    mutate(replic = recode(replic, "A" = "1", "B" = "2", "C" = "3")) |>
    unite("sample_id", c(donor_id, replic), sep = "_")

tmp <- 
    inner_join(distinct(fdr01, variant_id), ase_res) |>
    mutate(total = ref_count + alt_count) |>
    filter(fdr <= 0.01, total >= 40) |>
    distinct(variant_id)

ase_plot_data <- 
    inner_join(distinct(fdr01, variant_id), ase_res) |>
    inner_join(tmp) |>
    pivot_longer(ref_count:alt_count, names_to = "allele", values_to = "counts") |>
    mutate(allele = str_remove(allele, "_count")) |>
    mutate(sample_id = factor(sample_id, levels = sort(unique(sample_id)))) |>
    arrange(sample_id, stim, variant_id)

ase_plot_data_i <-
    ase_plot_data |>
    filter(variant_id == unique(variant_id)[1])

p <- 
    ggplot(ase_plot_data_i, 
	   aes(x = counts, y = allele)) +
    geom_col(aes(fill = stim)) +
    scale_x_continuous(breaks = scales::pretty_breaks(2)) +
    scale_fill_manual(values = stim_colors) +
    facet_grid2(sample_id~stim, scales = "free_x", independent = "x", drop = FALSE) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 11),
	  strip.text.x = element_text(size = 11),
	  strip.text.y = element_text(angle = 0),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(y = NULL) +
    guides(fill = "none")

ggsave("./plots/ase_time.png", p, height = 10, width = 5)

ase_plot_data_i |> 
    distinct(variant_id) |>
    inner_join(fdr01) |>
    arrange(fdr) |>
    mutate(fdr = fdr * 100,
	   fdr = round(fdr, 2)) |>
    print(n = Inf)




# Explore variants that do not replicate well
res_norand_df
ase_res

test_non <- 
    plot_data |> 
    filter((beta_stim.x > 0 & beta_stim.y < 0) |
	   (beta_stim.x < 0 & beta_stim.y > 0)) |>
    mutate(d = abs(beta_stim.x - beta_stim.y)) |>
    select(donor_id, stim, variant_id, replic.x, replic.y, d) |>
    mutate_at(vars(replic.x, replic.y), ~recode(., "A" = 1, "B" = 2, "C" = 3)) |> 
    mutate_at(vars(replic.x, replic.y), ~paste(donor_id, ., sep = "_")) |>
    select(donor_id, sample_1 = replic.x, sample_2 = replic.y, stim, variant_id, d) |>
    arrange(desc(d))


plot_data |>
    filter(stim == "TLR7", donor_id == "10059706", 
	   replic.x == "A", replic.y == "B",
	   beta_stim.x < 0, beta_stim.y > 0)


# Variant 1
var_x <- "chr7:74405694:C:G"
test_non |> filter(variant_id == var_x)

test_ase_data1 <-
    ase_res |>
    filter(sample_id %in% c("10059706_1", "10059706_2"), 
	   variant_id == var_x,,
	   stim %in% c("Day 0", "TLR7")) |>
    pivot_longer(ref_count:alt_count, names_to = "allele", values_to = "counts") |>
    mutate(allele = str_remove(allele, "_count"),
	   allele = toupper(allele)) |>
    select(sample_id, stim, variant_id, allele, counts, p_value)

test_ase_plot1 <-
    ggplot(test_ase_data1) +
    geom_col(aes(x = counts, y = allele, fill = stim)) +
    scale_fill_manual(values = stim_colors) +
    facet_grid(sample_id~stim) +
    theme_minimal() +
    theme(panel.grid.minor.x = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  plot.background = element_rect(fill = "white", color = "white")) +
    guides(fill = "none")

ggsave("./plots/test_ase_plot1.png", test_ase_plot1, height = 2, width = 5)

# Variant 2
ix <- 1

test_non |> slice(ix)

test_ase_data2 <- ase_res |>
    filter(sample_id %in% c("10059706_1", "10059706_2"),
	   variant_id == test_non$variant_id[ix]) |>
    pivot_longer(ref_count:alt_count, names_to = "allele", values_to = "counts") |>
    mutate(allele = str_remove(allele, "_count"),
	   allele = toupper(allele)) |>
    select(sample_id, stim, variant_id, allele, counts, p_value)

test_ase_plot2 <-
    ggplot(test_ase_data2) +
    geom_col(aes(x = counts, y = allele, fill = stim)) +
    scale_fill_manual(values = stim_colors) +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    facet_grid(sample_id~stim, scales = "free_x") +
    theme_minimal() +
    theme(panel.grid.minor.x = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  plot.background = element_rect(fill = "white", color = "white")) +
    guides(fill = "none")

ggsave("./plots/test_ase_plot2.png", test_ase_plot2, height = 2, width = 6)

# Variant 3
plot_data |>
    filter(donor_id == "10061340", 
	   stim == "BCR",
	   replic.x == "C", replic.y == "B",
	   beta_stim.x > 0, beta_stim.y < 0)

test_ase_data3 <- 
    ase_res |>
    filter(sample_id %in% c("10061340_2", "10061340_3"),
	   variant_id == "chr17:42700823:G:A", 
	   stim %in% c("Day 0", "BCR")) |>
    pivot_longer(ref_count:alt_count, names_to = "allele", values_to = "counts") |>
    mutate(allele = str_remove(allele, "_count"),
	   allele = toupper(allele)) |>
    select(sample_id, stim, variant_id, allele, counts, p_value)

test_ase_plot3 <-
    ggplot(test_ase_data3) +
    geom_col(aes(x = counts, y = allele, fill = stim)) +
    scale_fill_manual(values = stim_colors) +
    facet_grid(sample_id~stim) +
    theme_minimal() +
    theme(panel.grid.minor.x = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  plot.background = element_rect(fill = "white", color = "white")) +
    guides(fill = "none")

ggsave("./plots/test_ase_plot3.png", test_ase_plot3, height = 2, width = 5)



# Test: compare p_dev with p_reg
glm_df <- 
    "./results_glm/glm_res_df.tsv" |>
    read_tsv() |>
    filter(grepl("^Day 0", stim)) |>
    mutate(replic = recode(replic, "A" = "1", "B" = "2", "C" = "3"),
	   stim = str_remove(stim, "Day 0-")) |>
    arrange(donor_id, replic) |>
    unite("sample_id", c(donor_id, replic), sep = "_") |>
    mutate(sample_id = fct_inorder(sample_id),
	   stim = factor(stim, levels = c("BCR", "TLR7", "DN2")))
    
p_vals_plot <- 
    ggplot(glm_df, aes(x = p_dev, y = p_reg)) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(color = stim), alpha = .25, size = .5) +
    scale_x_continuous(limits = c(0, 1),
		       breaks = c(0, 1),
		       labels = c("0", "1")) +
    scale_y_continuous(limits = c(0, 1),
		       breaks = c(0, 1),
		       labels = c("0", "1")) +
    scale_color_manual(values = stim_colors) +
    facet_grid(sample_id~stim) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  plot.background = element_rect(color = "white", fill = "white")) +
    guides(color = "none") +
    labs(x = "P-value (Deviance)", y = "P-value (regression)")

ggsave("./plots/glm_pvals_comparison.png", p_vals_plot, height = 9, width = 3)
