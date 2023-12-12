library(tidyverse)
library(qvalue)
library(extrafont)
library(patchwork)
library(ggbeeswarm)


ase_data <- read_tsv("ase_data.tsv", col_types = "ffccccii")

#imb_df <- ase_data |>
#    mutate(imb = abs(0.5 - (ref_count/total_count))) |>
#    select(sample_id, stim, variant_id, total_count, imb)

ase_res <- 
    ase_data |>
    mutate(total = refCount + altCount) |>
    mutate(p_value = map2_dbl(refCount, total, 
			      ~binom.test(.x, .y, p = .5, alternative = "two.sided")$p.value),
	   q_value = qvalue(p_value)$qvalues)

# Plot parameters
stim_colors <- 
    c("Day 0" = "#898e9f",
      "BCR" = "#003967",
      "TLR7" = "#637b31",
      "DN2" = "#a82203")

text_size_big <- 9
text_size_small <- 8

# Plots
summary_1 <- 
    ase_res |>
    group_by(stim) |>
    summarise(n = n_distinct(var_id)) |>
    ungroup()


summary_1_plot <- 
    summary_1 |>
    mutate(labtext = scales::comma(n)) |>
    ggplot(aes(x = stim, y = n)) +
    geom_col(aes(fill = stim)) +
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
	 title = "Number of different variants tested\nper condition") +
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
    geom_col(aes(fill = stim), position = "dodge", color = "black") +
    scale_fill_manual(values = stim_colors) + 
    scale_y_continuous(limits = c(0, 8e3)) +
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
    geom_col(aes(fill = stim), position = "dodge", color = "black") +
    scale_fill_manual(values = stim_colors) + 
    scale_y_continuous(limits = c(0, 4e3)) +
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
    geom_quasirandom(aes(color = stim), method = "smiley", width = .2) +
    scale_y_continuous(limits = c(0, 30e3), label = scales::comma) +
    scale_color_manual(values = stim_colors) +
    theme_minimal() +
    theme(
	  axis.text = element_text(family = "Arial", size = text_size_small),
	  axis.title = element_text(family = "Arial", size = text_size_big),
	  plot.background = element_rect(color = "white", fill = "white"),
	  plot.title = element_text(family = "Arial", size = text_size_big),
	  plot.title.position = "plot"
	  ) +
    labs(x = NULL, y = "Number of variants",
	 title = "Total number of variants tested\nper sample") +
    guides(color = "none")

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
		     method = "smiley", width = .2, size = .75) +
    scale_y_continuous(limits = c(0, 1400), 
		       breaks = seq(0, 1400, 200),
		       label = scales::comma) +
    scale_color_manual(values = stim_colors) +
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
    guides(color = "none")


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
		     method = "smiley", width = .2, size = .75) +
    scale_y_continuous(limits = c(0, 800), 
		       breaks = seq(0, 800, 200)) +
    scale_color_manual(values = stim_colors) +
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
    guides(color = "none")



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

ggsave("./plots/qq.png", qq_plot, height = 12, width = 5, dpi = 600)

