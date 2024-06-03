library(tidyverse)
library(tidytext)
library(DESeq2)
library(extrafont)
library(furrr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(cowplot)

slice <- dplyr::slice
count <- dplyr::count
select <- dplyr::select

if (!file.exists("paper_plots")) dir.create("paper_plots")

stim_colors <- 
    "./figure_colors.txt" |>
    read_tsv(col_names = c("condition", "color")) |>
    mutate(condition = sub("_", " ", condition),
	   condition = paste0(condition, "hrs")) |>
    deframe()

dat <- read_rds("./bcell_lowinput/wgcna/data/gene_expression.rds")

sample_table <- 
    tibble(sample_id = colnames(dat$counts)) |>
    extract(sample_id, c("group"), "[^_]+_(.+)", remove = FALSE) |>
    column_to_rownames("sample_id")

dds <- 
    DESeqDataSetFromTximport(dat, sample_table, ~1) |>
    estimateSizeFactors() |>
    {function(x) x[, colSums(counts(x)) > 2e6]}() |>
    {function(x) x[rowSums(counts(x, normalized = TRUE) >= 5 ) >= 3, ]}() |>
    vst()

top_variable_genes <- 
    rowVars(assay(dds), useNames = TRUE) |> 
    enframe("gene_id", "var") |>
    arrange(desc(var)) |>
    slice(1:2000)

## Transpose expression matrix to use with WGCNA
count_matrix <- t(assay(dds)[top_variable_genes$gene_id, ])

# Run PCA
pca <- prcomp(count_matrix, center = TRUE, scale. = TRUE)

pc_scores <- as_tibble(pca$x, rownames = "sample_id")

pc_varexp <- 
    tibble(PC = paste0("PC", 1:length(pca$sdev)),
	   variance = pca$sdev^2) |>
    mutate(pct = round(variance/sum(variance) * 100, 1),
           pct = paste0("(", pct, "%)")) |>
    select(PC, pct) |>
    unite("lab", c(PC, pct), sep = " ", remove = FALSE)


# Plot
pca_df <- 
    pc_scores |>
    select(sample_id, PC1:PC2) |>
    separate(sample_id, c("donor_id", "stim", "timep"), sep = "_") |>
    mutate(condition = paste(stim, timep, sep = " "),
	   condition = factor(condition, levels = names(stim_colors)))

legend_df <- 
    stim_colors |>
    enframe("condition", "color") |>
    separate("condition", c("stim", "time"), sep = " ", remove = FALSE) |>
    filter(condition %in% pca_df$condition) |>
    mutate(time = str_remove(time, "hrs"),
	   time = fct_inorder(time),
	   stim = fct_inorder(stim),
	   stim = fct_rev(stim)) |>
    select(stim, time, condition)

legend_plot <-
    ggplot(legend_df, aes(x = time, y = stim, fill = condition)) +
    geom_point(size = 2.5, shape = 21, stroke = .25) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  axis.text.y = element_text(margin = margin(r = -.5, l = .25, unit = "lines")),
	  axis.title.y =  element_blank(),
	  panel.grid = element_blank(),
	  plot.title = element_text(size = 9),
	  plot.margin = unit(c(0, 0, 0, 0), "lines")) +
    labs(x = "hours") +
    guides(fill = "none")

pca_plot <- 
    ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(fill = condition), 
	       size = 2.5, shape = 21, stroke = .25) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  panel.grid = element_blank(),
	  plot.title = element_text(size = 9)) +
    guides(fill = "none") +
    labs(x = pc_varexp$lab[1], y = pc_varexp$lab[2])

pca_title <- 
    ggdraw() + 
    draw_label(
	       "PCA shows separation of stimuli and time points",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(0, 1.25, 0, 1.25, unit = "lines"))

pca_grid <- 
    plot_grid(
	      pca_title,
	      plot_grid(
			pca_plot, 
			plot_grid(NULL, legend_plot, NULL, ncol = 1, rel_heights = c(.1, 1, .8)), 
			nrow = 1, 
			rel_widths = c(1, .6)
			),
	      ncol = 1, 
	      rel_heights = c(.1, 1)
    )


# WGCNA
module_sizes <-
    list.files("./bcell_lowinput/wgcna/data",
	       pattern = "_network\\.rds",
	       full.names = TRUE) |>
    {function(x) setNames(x, str_extract(basename(x), "[^_]+"))}() |>
    map_dfr(~read_rds(.) |>
	    {function(x) table(x$colors)}() |>
	    enframe("module", "n") |>
	    filter(module != "grey") |>
	    mutate(n = as.integer(n)) |>
	    arrange(desc(n)) |>
	    mutate(ix = fct_inorder(module),
		   ix = as.integer(ix)),
	    .id = "stim")

eigengenes_df <-
    list.files("./bcell_lowinput/wgcna/data",
	       pattern = "eigen\\.tsv",
	       full.names = TRUE) |>
    {function(x) setNames(x, str_extract(basename(x), "[^_]+"))}() |>
    map_dfr(~read_tsv(.) |>
	    separate(sample_name, 
		     c("donor_id", "stim", "time"), 
		     sep = "_") |>
	    pivot_longer(starts_with("ME"), names_to = "module") |>
	    mutate(module = str_remove(module, "^ME"),
		   time = parse_number(time),
		   time = factor(time, levels = sort(unique(time)))) |>
	    filter(module != "grey"),
	    .id = "stim") |>
    left_join(module_sizes, join_by(stim, module)) |>
    select(stim, donor_id, module_ix = ix, module, time, value) |>
    arrange(stim, module_ix, donor_id, time)

stim_colors_single <-
    stim_colors |>
    {function(x) keep(x, !grepl("Unstim|IL4", names(x)))}() |>
    {function(x) keep(x, grepl("48hrs", names(x)))}() |>
    {function(x) setNames(x, str_extract(names(x), "[^ ]+"))}()

wgcna_plot <- 
    ggplot(data = eigengenes_df, 
       aes(x = time, y = value)) +
    geom_line(aes(group = donor_id),
	      linewidth = .5, alpha = .5, 
	      show.legend = FALSE) +
    facet_grid(factor(stim, levels = names(stim_colors_single)) ~ factor(module_ix),
	       scales = "free_y") +
    theme_minimal() +
    theme(text = element_text(size = 8),
	  axis.text = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  panel.grid.major.x = element_line(color = "grey90", linewidth = .2),
	  panel.grid.major.y = element_blank(),
	  strip.text.y = element_text(angle = 0, hjust = 0),
	  panel.spacing = unit(0.1, "lines")
	  ) +
    labs(x = "hours", y = "Eigengene expression")

wgcna_title <- 
    ggdraw() + 
    draw_label(
	       "Gene programs in stimulated B cells",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(0, 1, 0, 1, unit = "lines"))

wgcna_grid <- 
    plot_grid(
	      wgcna_title, wgcna_plot,
	      ncol = 1, rel_heights = c(.1, 1)
    )

fig2_top_grid <-
    plot_grid(pca_grid, NULL, wgcna_grid, 
	      nrow = 1, rel_widths = c(1, 0.05, .85),
	      labels = c("A)", "B)"), label_size = 10, label_x = -0.0125)


# Fig 2 bottom

dn2_modules_df <-
    eigengenes_df |>
    filter(stim == "DN2") |>
    group_by(stim, module = module_ix, time) |>
    summarise(value = mean(value)) |>
    ungroup() |>
    mutate(module = fct_inorder(as.character(module)))

dn2_modules_plot <- 
    ggplot() +
    geom_line(data = dn2_modules_df,
	      aes(x = time, y = value, group = module),
	      linewidth = 1.25) +
    facet_wrap(~module, scale = "free_y", ncol = 1) +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  axis.text.y = element_blank(),
	  axis.ticks.y = element_blank(),
	  panel.grid = element_blank(),
	  strip.text = element_text(size = 9,
				    margin = margin(t = 0, b = 0)),
	  strip.background = element_rect(fill = "grey80", color = NA)) +
    labs(x = "hours", y = "Average Eigengene expression")


kim_df <- 
    "./bcell_lowinput/wgcna/data/DN2_kim.tsv" |>
    read_tsv() |>
    group_by(module) |>
    top_n(5, kim) |>
    ungroup() |>
    mutate(stim = "DN2") |>
    left_join(module_sizes, join_by(stim, module)) |>
    select(module = ix, gene_name, kim) |>
    arrange(module, desc(kim))

kim_plot <-
    ggplot(data = kim_df,
	   aes(x = kim, 
	       y = reorder_within(gene_name, by = kim, within = module))) +
    geom_col() +
    scale_x_continuous(breaks = c(.5, 1),
		       labels = c("0.5", "1")) +
    scale_y_reordered(position = "right") +
    facet_wrap(~module, scale = "free_y", ncol = 1) +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  axis.text.y.right = element_text(hjust = 0, 
					   margin = margin(l = -.5, unit = "lines")),
	  panel.grid = element_blank(),
	  strip.text = element_blank()) +
    labs(x = "kIM", y = NULL) + 
    coord_cartesian(xlim = c(0.5, 1))

# GO
go_res <-
    "./bcell_lowinput/wgcna/data/DN2_go.tsv" |>
    read_tsv() |>
    group_by(module) |>
    top_n(5, -log10(pvalue)) |>
    ungroup() |>
    mutate(stim = "DN2",
	   Description = str_trunc(Description, width = 30),
	   gene_r = map_dbl(GeneRatio, ~eval(parse(text = .)))) |>
    left_join(module_sizes, join_by(stim, module)) |>
    select(module = ix, Description, pvalue, gene_r) |>
    mutate(module = factor(module, levels = unique(kim_df$module))) 

go_plot <- 
    ggplot(data = go_res, 
       aes(x = "1", 
	   y = reorder_within(Description, by = -log10(pvalue), within = module))) +
    geom_point(aes(fill = -log10(pvalue), size = gene_r), 
	       shape = 21, stroke = .5) +
    scale_y_reordered(position = "right") +
    scale_fill_gradient(low = "beige", high = "firebrick") +
    scale_size(range = c(.5, 3.5)) +
    facet_wrap(~module, drop = FALSE,
	       ncol = 1, scale = "free_y") +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  axis.text.x = element_blank(),
	  legend.text = element_text(size = 9, margin = margin(l = -0, unit = "lines")),
	  legend.title = element_text(size = 9),
	  strip.text = element_blank(),
	  panel.grid = element_blank()) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colorbar("logP:", barheight = 5, barwidth = .5),
	   size = guide_legend("Gene\nRatio:")) +
    coord_cartesian(clip = "off")

dn2_title <- 
    ggdraw() + 
    draw_label(
	       "DN2 condition: gene programs are enriched\nin different biological processes",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(-1, 2, -1, 2, unit = "lines"))

dn2_grid <-
    plot_grid(
	      dn2_title,
	      plot_grid(dn2_modules_plot, kim_plot, NULL, go_plot,
			nrow = 1, align = "h", rel_widths = c(.6, .5, .1, 1)),
	      ncol = 1, rel_heights = c(.1, 1))

fig2_bottom_grid <-
    plot_grid(dn2_grid, NULL, 
	      nrow = 1, rel_widths = c(1, .6),
	      labels = c("C)", "D)"), label_size = 10, label_x = -0.0125)

fig2_grid <- 
    plot_grid(fig2_top_grid, fig2_bottom_grid, ncol = 1, rel_heights = c(0.4, 1)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./paper_plots/fig2.png", fig2_grid, width = 6.5, height = 2.25 + 5.5)

