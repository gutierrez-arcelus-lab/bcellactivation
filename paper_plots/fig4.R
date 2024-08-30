library(tidyverse)
library(ggrepel)
library(scico)
library(cowplot)
library(glue)
library(tidytext)

# colors
stim_colors <- 
    read_tsv("../figure_colors.txt", col_names = c("stim", "color")) |>
    mutate(stim = str_replace(stim, "_", " "),
	   stim = paste0(stim, "h")) |>
    deframe()

# meta data
donor_ids <- 
    "../atacseq/samplesheet.csv" |>
    read_csv() |>
    mutate(donor_id = basename(fastq_1)) |>
    extract(donor_id, "donor_id", "[^_]+_([^_]+)_.+") |>
    mutate(replicate = paste0("REP", replicate),
	   stim = str_replace(sample, "_", " "),
	   stim = str_replace(stim, "unst", "Unstim"),
	   stim = paste0(stim, "h")) |>
    select(stim, donor_id, donor = replicate) |>
    distinct() |>
    filter(donor_id != "3donors")

# Figure A ####################################################################
pca_file <- "../atacseq/results_deseq2/pcadata_5000peaks.rds"

pca_5000_df <- 
    read_rds(pca_file) |>
    mutate(stim = str_replace(condition, "_", " "),
	   stim = str_replace(stim, "unst", "Unstim"),
	   stim = paste0(stim, "h")) |>
    extract(name, "donor", ".+(REP\\d)") |>
    select(stim, donor, PC1, PC2) |>
    left_join(donor_ids) |>
    mutate(stim = factor(stim, levels = names(stim_colors)))

perc_var <-
    read_rds(pca_file) |>
    attr("percentVar") |>
    {function(x) round(x * 100)}()

fig_a <- 
    ggplot(pca_5000_df, aes(PC1, PC2)) +
    geom_point(aes(fill = stim), shape = 21, size = 3.5, stroke = .5) +
    scale_fill_manual(values = stim_colors) +
    theme_bw() +
    theme(
	  axis.text = element_text(size = 8),
	  axis.title = element_text(size = 9),
	  legend.text = element_text(size = 8, margin = margin(l = 0)),
	  legend.title = element_text(size = 9),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  panel.grid.major.x = element_line(color = "grey96"),
	  panel.grid.major.y = element_line(color = "grey96"),
	  #legend.justification.top = "left",
	  legend.key.spacing.y = unit(.5, "pt"),
	  legend.box.spacing = unit(1, "pt")
	  ) +
    labs(x = sprintf("PC1: %s%% variance", perc_var[1]),
	 y = sprintf("PC2: %s%% variance", perc_var[2]),
	 fill = "Stim:")

fig_a_title <- 
    ggdraw() + 
    draw_label(
	       "PCA shows separation of conditions",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_a_grid <- plot_grid(fig_a_title, fig_a, ncol = 1, rel_heights = c(.1, 1))




# Figure B ####################################################################
da_files <- 
    list.files("../atacseq/results_deseq2",
	       pattern = ".+vs.+\\.tsv$",
	       full.names = TRUE)

da_files <- setNames(da_files, basename(da_files) |> str_remove("\\.tsv"))

da_data <-
    map_dfr(da_files, read_tsv, .id = "comparison") |>
    filter(!is.na(padj)) |>
    separate(comparison, c("stim1", "stim2"), sep = "vs")

da_positive <- 
    da_data |>
    filter(log2FoldChange > 0, padj <= 0.01)

da_negative <- 
    da_data |>
    filter(log2FoldChange < 0, padj <= 0.01)

da_summ <- 
    bind_rows("positive" = da_positive, "negative" = da_negative, .id = "comparison") |>
    mutate(comparison = fct_inorder(comparison)) |>
    count(comparison, stim1, stim2) |>
    mutate_at(vars(stim1, stim2), 
	      ~str_replace(., "_", " ") |> 
	      str_replace("unst", "Unstim") |>
	      paste0("h") |> 
	      factor(levels = names(stim_colors))) |>
   arrange(comparison, stim2, stim1) |>
   mutate(lab = n,
	  n = case_when(comparison == "positive" ~ n,
			comparison == "negative" ~ n * -1),
	  stim_a = case_when(comparison == "positive" ~ stim1,
			     comparison == "negative" ~ stim2),
	  stim_b = case_when(comparison == "positive" ~ stim2,
			     comparison == "negative" ~ stim1)) |>
    select(comparison, stim_a, stim_b, n, lab)

da_plot <- 
    ggplot(da_summ, aes(x = stim_b, y = stim_a)) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = scales::comma(lab)), 
	      color = "white", fontface = "bold", size = 7, size.unit = "pt") +
    scale_x_discrete(position = "top") +
    scale_fill_scico(palette = "vik", na.value = "white") +
    theme_minimal() + 
    theme(
	  panel.grid = element_blank(),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.ticks.y = element_blank(),
	  legend.position = "none",
	  plot.margin = margin(0, 0, 0, 0, "cm")
    )

cond_vertical <- 
    ggplot(da_summ |> distinct(stim_a),
	   aes(x = factor(1), y = stim_a)) +
    geom_tile(aes(fill = stim_a)) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.ticks.y = element_blank(),
	  plot.margin = margin(0, 0, 0, 0)) +
    guides(fill = "none") +
    coord_cartesian(xlim = c(1, 1))

cond_horizontal <- 
    ggplot(da_summ |> distinct(stim_b),
	   aes(y = factor(1), x = stim_b)) +
    geom_tile(aes(fill = stim_b)) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.ticks.y = element_blank(),
	  plot.margin = margin(0, 0, 0, 0, "lines")) +
    guides(fill = "none") +
    coord_cartesian(ylim = c(1, 1))


fig_b <- 
    plot_grid(NULL, cond_horizontal, cond_vertical, da_plot,
	      ncol = 2, nrow = 2, align = "hv", rel_heights = c(.085, 1), rel_widths = c(0.05, 1))


fig_b_title <- 
    ggdraw() + 
    draw_label(
	       "Number of differentially accessible peaks",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_b_grid <- 
    plot_grid(fig_b_title, fig_b, ncol = 1, rel_heights = c(.1, 1)) +
    theme(plot.margin = margin(r = 5.5 * 2, unit = "pt"))
	


# Figure C ####################################################################
stims <- read_lines("../atacseq/homer/stims.txt")

homer_df <- 
    glue("../atacseq/homer/results/{stims}/knownResults.txt") |>
    setNames(stims) |>
    map_dfr(~read_tsv(.) |> 
	    janitor::clean_names() |>
	    {function(x) setNames(x, str_remove(names(x), "_of_\\d+$"))}(),
	    .id = "stim") |>
    mutate_at(vars(starts_with("percent_of_")), parse_number) |>
    mutate(stim = factor(stim, levels = stims),
	   fc = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) |>
    filter(!grepl("bias", motif_name, ignore.case = TRUE)) |>
    mutate(motif_name = str_replace(motif_name, "NR/Ini-like", "NR;Ini-like")) |>
    separate(motif_name, c("motif", "dataset", "db"), sep = "/") |>
    separate(motif, c("motif", "tf_family"), sep = "\\(") |>
    mutate(tf_family = str_remove(tf_family, "\\)$")) |>
    mutate(log10p = log_p_value/log(10)) |>
    select(stim, motif, tf_family, dataset, 
	   pct_target =  percent_of_target_sequences_with_motif,
	   pct_bg = percent_of_background_sequences_with_motif,
	   fc, log10p, q_value = q_value_benjamini)

homer_top25 <- 
    homer_df |>
    filter(q_value <= 0.05) |>
    group_by(stim) |>
    top_n(25, -log10p) |>
    ungroup()


homer_top25_inall <- 
    homer_df |>
    inner_join(distinct(homer_top25, motif, tf_family)) |>
    arrange(log10p) |>
    mutate(motif = fct_inorder(motif))

homer_heat <- 
    ggplot(homer_top25_inall, 
	   aes(x = motif, y = stim)) +
    geom_tile(aes(fill = log2(fc))) +
    geom_text(data = filter(homer_top25_inall, q_value <= 0.01), 
	      aes(x = motif, y = stim, label = "*"), 
	      size = 7, size.unit = "pt", nudge_y = -.1) +
    scale_y_discrete(position = "left") +
    scale_fill_scico(palette = "vik", midpoint = 0) +
    facet_grid(cols = vars(tf_family), scales = "free", space = "free", switch = "x") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1, 
				     margin = margin(0, 0, 0, 0)),
	  axis.text.y = element_text(size = 8),
	  axis.title = element_blank(),
	  strip.text.x.bottom = element_blank(),
	  panel.spacing.x = unit(0.1, "pt"),
	  legend.position = "right",
	  legend.text = element_text(size = 8),
	  legend.title = element_text(size = 8),
	  legend.box.spacing = unit(-.25, "lines"),
	  plot.margin = margin(0, 0, 0, 0),
	  plot.background = element_rect(fill = "white", color = "white")) +
    guides(fill = guide_colorbar(barwidth = .25, barheight = 1.5)) +
    labs(fill = NULL)



homer_plot <- 
    ggplot(homer_top25, 
	   aes(x = -log10p, 
	       y = reorder_within(motif, by = -log10p, within = stim))) +
    geom_point(aes(fill = log2(fc)), shape = 21, size = 2, stroke = .2) +
    scale_x_continuous(limits = function(x) c(0, max(x)),
		       expand = expansion(mult = .1),
		       breaks = function(x) c(0, max(x)),
		       labels = function(x) round(x)) +
    scale_y_reordered() +
    scale_fill_gradient(low = "white", high = "black",
			limits = function(x) c(0, max(x)),
			breaks = function(x) seq(0, max(x), length.out = 3),
			labels = function(x) round(x)) +
    facet_wrap(~stim, nrow = 1, scales = "free") +
    theme_minimal() +
    theme(axis.text = element_text(size = 8),
	  axis.title.x = element_text(size = 9, margin = margin(t = -.5, unit = "lines")),
	  axis.title.y = element_blank(),
	  strip.text = element_text(size = 9),
	  strip.clip = "off",
	  legend.box.spacing = unit(.1, "lines"),
	  legend.text = element_text(size = 8),
	  legend.title = element_text(size = 8),
	  panel.grid.minor = element_blank(),
	  plot.margin = margin(0, 5.5, 0, 5.5, "pt")
	  ) +
    guides(fill = guide_colorbar(barwidth = .25))


fig_c <- plot_grid(homer_plot, NULL, homer_heat, ncol = 1, rel_heights = c(1, .05, .4))

fig_c_title <- 
    ggdraw() + 
    draw_label(
	       "Motif enrichment of DA peaks in a condition in respect to Unstimulated 24 hours",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_c_grid <- plot_grid(fig_c_title, fig_c, ncol = 1, rel_heights = c(.1, 1),
			labels = "c", label_size = 12)

# Figure D ####################################################################

# Figure E ####################################################################


# Final panel
top_grid <- 
    plot_grid(fig_a_grid, NULL, fig_b_grid, nrow = 1, rel_widths = c(1, .05, .85),
	      labels = c("a", "", "b"), label_size = 12)

final_grid <- 
    plot_grid(top_grid, NULL, fig_c_grid, ncol = 1, rel_heights = c(.5, .01, 1)) +
    theme(panel.background = element_rect(color = "white", fill = "white"))


ggsave("fig4.png", final_grid, width = 6.5, height = 6.5)

