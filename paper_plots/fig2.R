library(tidyverse)
library(tidytext)
library(DESeq2)
library(extrafont)
library(furrr)
library(cowplot)
library(ggbeeswarm)
library(glue)
library(ggh4x)

slice <- dplyr::slice
count <- dplyr::count
select <- dplyr::select
filter <- dplyr::filter

stim_colors <- 
    "../figure_colors.txt" |>
    read_tsv(col_names = c("condition", "color")) |>
    mutate(condition = sub("_", " ", condition),
	   condition = paste0(condition, "hrs")) |>
    deframe()

stim_colors_0_names <- 
    keep(names(stim_colors), !grepl("Unstim", names(stim_colors))) |>
    str_remove(" \\d+hrs$") |>
    unique() |>
    paste("0hrs")

stim_colors_0 <- rep(stim_colors[[1]], length(stim_colors_0_names))
names(stim_colors_0) <- stim_colors_0_names

stim_colors <- c(stim_colors, stim_colors_0)

# Fig A PCA
dat <- read_rds("../bcell_lowinput/wgcna/data/gene_expression.rds")

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
    mutate(timep = factor(timep, levels = paste0(c(0, 4, 24, 48, 72), "hrs")),
	   condition = paste(stim, timep, sep = " "),
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
    select(stim, time, condition) |>
    mutate(stim = paste0(stim, ":"),
	   stim = fct_inorder(stim))

legend_plot <-
    ggplot(legend_df, aes(x = time, y = 1)) +
    geom_point(aes(fill = condition, shape = time), 
	       size = 2.5, stroke = .25) +
    geom_point(data = filter(legend_df, time == 0),
	       size = 1.5, shape = 4) +
    scale_fill_manual(values = stim_colors) +
    scale_shape_manual(values = c(21, 21, 23, 22, 24)) +
    facet_wrap(~stim, ncol = 1) +
    theme_minimal() +
    theme(text = element_text(size = 8),
	  axis.text.y = element_blank(),
	  axis.title.y =  element_blank(),
	  axis.ticks.x = element_line(linewidth = .25),
	  axis.ticks.y = element_blank(),
	  panel.grid = element_blank(),
	  panel.border = element_rect(color = NA, fill = NA),
	  panel.spacing = unit(0, "lines"),
	  strip.clip = "off",
	  strip.text = element_text(margin = margin(t = 0, b = 0)),
	  strip.background = element_rect(color = NA, fill = "grey90"),
	  plot.margin = margin(0, 2, 0, -1, unit = "pt")
	  ) +
    labs(x = "hours") +
    guides(fill = "none", shape = "none")

pca_plot <- 
    ggplot(pca_df, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, color = "grey94") +
    geom_vline(xintercept = 0, color = "grey94") +
    geom_point(data = filter(pca_df, timep == "0hrs"), 
	       aes(fill = condition, shape = timep), 
	       size = 2, stroke = .15) +
    geom_point(data = filter(pca_df, timep == "0hrs"),
	       size = 1.5, shape = 4) +
    geom_point(data = filter(pca_df, timep != "0hrs"), 
	       aes(fill = condition, shape = timep), 
	       size = 2, stroke = .15) +
    scale_fill_manual(values = stim_colors) +
    scale_shape_manual(values = c(21, 21, 23, 22, 24)) +
    theme_minimal() +
    theme(
	  axis.text = element_blank(),
	  axis.title = element_text(size = 9),
	  panel.grid = element_blank(),
	  plot.margin = margin(0, 0, 0, .5, unit = "lines"),
	  plot.title = element_text(size = 9)) +
    guides(fill = "none", shape = "none") +
    labs(x = pc_varexp$lab[1], y = pc_varexp$lab[2])

pca_inset <- 
    ggplot(pca_df |> mutate(i = stim %in% c("BCR", "BCR-TLR7", "DN2")),
	   aes(PC1, PC2)) +
    geom_point(aes(color = i), 
	       size = .5) +
    scale_color_manual(values = c("TRUE" = "midnightblue", "FALSE" = "grey90")) +
    theme_bw() +
    theme(axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.ticks = element_blank(),
	  panel.grid = element_blank(),
	  panel.border = element_rect(color = "grey"),
	  panel.background = element_rect(fill = "transparent"),
	  plot.background = element_rect(fill = "transparent", color = NA),
	  plot.title = element_text(size = 7, margin = margin(b = 0)),
	  plot.margin = margin(0, 0, 0, 0)) +
    guides(color = "none") +
    labs(title = "Includes BCR\nstimulation:")

pca_grid <- 
    ggdraw(pca_plot) +
    draw_plot(pca_inset, 0.15, 0.1, .3, .35)

fig_a_title <- 
    ggdraw() + 
    draw_label(
	       "PCA shows separation of stimuli and\ntime points",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 2.5, unit = "lines"))

fig_a_grid <- 
    plot_grid(
	      fig_a_title,
	      plot_grid(
			pca_plot, 
			plot_grid(legend_plot, NULL, pca_inset, ncol = 1, rel_heights = c(1, .05, .5)), 
			nrow = 1, 
			rel_widths = c(1, .3)
			),
	      ncol = 1, 
	      rel_heights = c(.1, 1),
	      labels = c("a", NULL), label_size = 12
    )



# Fig B #######################################################################
all_stims <-
    names(stim_colors) |>
    str_split(" ") |>
    map_chr(1) |>
    unique()

diff_expr <- 
    "../bcell_lowinput/results/edger/diff_expr_all_times.tsv" |>
    read_tsv()

diff_expr_summ <-
    diff_expr |>
    group_by(group1, group2) |>
    summarise(n = sum(!is.na(gene_id))) |>
    ungroup() |>
    separate(group1, c("stim1", "t1"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    separate(group2, c("stim2", "t2"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    mutate_at(vars(stim1, stim2), ~recode(., "BCR_TLR7" = "BCR-TLR7")) |>
    mutate_at(vars(stim1, stim2), ~factor(., levels = all_stims)) |>
    complete(group1, group2, fill = list(n = NA)) |>
    arrange(stim1, stim2, t1, t2) |>
    mutate_at(vars(group1, group2), ~factor(., levels = unique(c(group1, group2))))

segments_x <- 
    diff_expr_summ |>
    distinct(group1) |>
    rowid_to_column() |>
    separate(group1, c("stim", "time"), sep = "\\.") |>
    select(-time) |>
    group_by(stim) |>
    slice_max(rowid) |>
    ungroup() |>
    filter(rowid != max(rowid)) |>
    mutate(rowid = rowid + .5) |>
    arrange(rowid)

segments_y <- 
    diff_expr_summ |>
    distinct(group2) |>
    rowid_to_column() |>
    separate(group2, c("stim", "time"), sep = "\\.") |>
    select(-time) |>
    group_by(stim) |>
    slice_max(rowid) |>
    ungroup() |>
    filter(rowid != max(rowid)) |>
    mutate(rowid = rowid + .5) |>
    arrange(rowid)

diff_plot <- 
    ggplot(data = diff_expr_summ, 
       aes(x = group1, y = group2)) +
    geom_tile(aes(fill = n), alpha = .85) +
    geom_vline(xintercept = segments_x$rowid, color = "white", linewidth = .5) +
    geom_hline(yintercept = segments_y$rowid, color = "white", linewidth = .5) +
    scale_fill_viridis_c(option = "inferno",
			 na.value = "white",
			 labels = scales::comma) +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  panel.grid.major = element_blank(),
	  panel.background = element_rect(fill = "transparent", color = NA),
	  plot.margin = margin(0, 0, 0, 0),
	  legend.position.inside = c(.9, .4)
	  ) +
    labs(fill = "DE\ngenes:") +
    guides(fill = guide_colorbar(barwidth = .5, barheight = 7, position = "inside")) +
    coord_cartesian(xlim = c(1, 28), ylim = c(1, 28))


b_axis_x_df <- 
    diff_expr_summ |>
    filter(!is.na(stim1)) |>
    distinct(group1, stim1, t1) |>
    mutate(condition1 = glue("{stim1} {t1}hrs")) |>
    select(group1, condition1)

b_axis_y_df <- 
    diff_expr_summ |>
    filter(!is.na(stim2)) |>
    distinct(group2, stim2, t2) |>
    mutate(condition2 = glue("{stim2} {t2}hrs")) |>
    select(group2, condition2)

b_axis_x_plot <-
    ggplot(data = b_axis_x_df, 
	   aes(y = factor(1), x = group1)) +
    geom_tile(aes(fill = condition1)) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(axis.text = element_blank(),
	  axis.title = element_blank(),
	  panel.grid = element_blank(), 
	  plot.margin = margin(0, 0, 0, 0)
	  ) +
    guides(fill = "none") +
    coord_cartesian(xlim = c(1, 28), ylim = c(1, 1))

b_axis_y_plot <-
    ggplot(data = b_axis_y_df, 
	   aes(x = factor(1), y = group2)) +
    geom_tile(aes(fill = condition2)) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(axis.text = element_blank(),
	  axis.title = element_blank(),
	  panel.grid = element_blank(), 
	  plot.margin = margin(0, 0, 0, 0)
	  ) +
    guides(fill = "none") +
    coord_cartesian(xlim = c(1, 1), ylim = c(1, 28))

fig_b_title <- 
    ggdraw() + 
    draw_label(
	       "Number of differentially expressed genes\nacross all conditions and time points",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))
		 
b_tmp <- 
    plot_grid(NULL, b_axis_x_plot, b_axis_y_plot, diff_plot,
	      ncol = 2, nrow = 2, align = "v", rel_heights = c(.05, 1), rel_widths = c(0.05, 1))

fig_b_grid <- 
    plot_grid(
	      fig_b_title,
	      b_tmp,
	      ncol = 1, 
	      rel_heights = c(.1, 1),
	      labels = c("b"), label_size = 12
    )

top_grid <- 
    plot_grid(fig_a_grid, NULL, fig_b_grid, nrow = 1, rel_widths = c(1, .1, 1))


# Fig C #######################################################################
edger_res <- 
    "../bcell_lowinput/results/edger/results.tsv" |>
    read_tsv()

cpm_df <- 
    "../bcell_lowinput/results/edger/cpm.tsv" |>
    read_tsv() |>
    group_by(sample_id, stim) |>
    nest() |>
    ungroup() |>
    separate(sample_id, c("sample_id", "dummy", "time"), sep = "_") |>
    mutate(hours = parse_number(time),
	   hours = factor(hours, levels = sort(unique(hours))),
           condition = paste(stim, time, sep = " "),
	   stim = factor(stim, levels = all_stims)) |>
    unnest(cols = data) |>
    select(sample_id, condition, stim, hours, gene_id, gene_name, obs_cpm, obs_logcpm)

pathway_colors <-
    keep(stim_colors, grepl("48hrs", names(stim_colors)))

names(pathway_colors) <- str_remove(names(pathway_colors), " 48hrs")

pathway_colors <- c("IL4" = "goldenrod4", pathway_colors)

plot_timecourse <- function(counts_df) {

    ggplot(data = counts_df, aes(x = hours, y = obs_cpm)) +
    geom_quasirandom(aes(fill = condition),
		     method = "smiley", width = .2, 
		     shape = 21, stroke = .1, size = 1.25, alpha = .4) +
    geom_line(aes(group = stim, color = stim), 
	      stat = "smooth", method = "loess", span = 1, se = FALSE, 
	      linewidth = .7) +
    scale_y_continuous(expand = expansion(mult = c(0.25, 0.05)),
		       breaks = c(0, ceiling(max(counts_df$obs_cpm)))) +
    scale_color_manual(values = pathway_colors) + 
    scale_fill_manual(values = stim_colors) + 
    facet_grid(gene_name~stim, scale = "free", space = "free_x") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.spacing.x = unit(0, "null"),
	  strip.text = element_blank(),
	  strip.clip = "off",
	  axis.text.x = element_blank(),
	  axis.text.y = element_text(size = 8),
	  plot.title = element_text(size = 8, hjust = .5, 
				    margin = margin(b = 0)),
	  legend.position = "none",
	  plot.margin = margin(2, 2, 0, 2, "pt")) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL, title = unique(counts_df$gene_name))
}

fig_c_title <- 
    ggdraw() + 
    draw_label(
	       "B cell activation genes show distinct patterns of expression across conditions",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(b = .5, l = 4, unit = "lines"))
       
fig_c_grid <- 
    plot_grid(NULL,
	      plot_grid(plot_timecourse(filter(cpm_df, gene_name == "AICDA")),
			plot_timecourse(filter(cpm_df, gene_name == "FCER2")),
			plot_timecourse(filter(cpm_df, gene_name == "CD69")),
			ncol = 1),
	      NULL,
	      plot_grid(plot_timecourse(filter(cpm_df, gene_name == "CD86")),
			plot_timecourse(filter(cpm_df, gene_name == "TBX21")),
			plot_timecourse(filter(cpm_df, gene_name == "SLAMF7")),
			ncol = 1),
	      NULL,
	      nrow = 1, rel_widths = c(0.05, 1, 0.05, 1, 0.05)) +
    draw_label("Counts per million", x = -.001, y = 0.5, vjust = 1.5, angle = 90, size = 9)
       
fig_c <- 
    plot_grid(fig_c_title, fig_c_grid, ncol = 1, rel_heights = c(.1, 1),
	      labels = "c", label_size = 12, label_y = 1.75)

# Fig D #######################################################################
module_sizes <-
    read_tsv("../bcell_lowinput/wgcna/data/DN2_modules.tsv") |>
    count(module) |>
    arrange(desc(n)) |>
    filter(module != "grey") |>
    rowid_to_column("ix") |>
    select(module, ix) |>
    mutate(ix = paste("Module", ix)) |>
    mutate(ix = fct_inorder(ix))

module_colors <- setNames(module_sizes$module, module_sizes$module)
module_colors["yellow"] <- "goldenrod2"
module_colors["green"] <- "forestgreen"

eigengenes_df <-
    "../bcell_lowinput/wgcna/data/DN2_eigen.tsv" |>
    read_tsv() |>
    select(-MEgrey) |>
    separate(sample_name, 
	     c("donor_id", "stim", "time"), 
	     sep = "_") |>
    pivot_longer(starts_with("ME"), names_to = "module") |>
    mutate(module = str_remove(module, "^ME"),
	   time = parse_number(time),
	   time = factor(time, levels = sort(unique(time)))) |>
    left_join(module_sizes, join_by(module)) |>
    select(donor_id, module_ix = ix, module, time, value) |>
    arrange(module_ix, donor_id, time)

modules_df <-
    eigengenes_df |>
    group_by(module, module_ix, time) |>
    summarise(value = mean(value)) |>
    ungroup() |>
    arrange(module_ix, time)

#kim_df <- 
#    "../bcell_lowinput/wgcna/data/DN2_kim.tsv" |>
#    read_tsv() |>
#    left_join(filter(module_sizes, stim == "DN2"), join_by(module)) |> 
#    select(gene_id, gene_name, module, module_ix = ix, kim) |>
#    arrange(module_ix, desc(kim)) |>
#    group_by(module) |>
#    top_n(5, kim) |>
#    ungroup()

fig_d <-
    ggplot(modules_df) +
    geom_line(aes(x = time, y = value, group = module, color = module),
	      linewidth = 1.25) +
    scale_color_manual(values = module_colors) +
    facet_wrap(~module_ix, ncol = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8),
	  axis.text.y = element_text(size = 8),
	  axis.title.x = element_blank(),
	  axis.title.y = element_text(size = 9),
	  panel.grid.major.x = element_line(linetype = 2, color = "black", linewidth = .1),
	  panel.grid.major.y = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  strip.text = element_text(size = 8.5, margin = margin(b = 0)),
	  strip.clip = "off",
	  plot.margin = margin(.5, 0, 0, 0, unit = "lines")) +
    guides(color = "none") +
    labs(y = "Eigengene expression")


#kim_plot <-
#    ggplot(data = kim_df |> 
#	   arrange(module_ix, kim) |>
#	   mutate(gene_name = fct_inorder(gene_name)),
#	   aes(x = " ", y = gene_name)) +
#    geom_point(aes(fill = kim), shape = 21, size = 3, show.legend = FALSE) +
#    scale_x_discrete(drop = FALSE, expand = c(0, 0), limits = c(" ", " ")) +
#    scale_y_discrete(position = "right") +
#    scale_fill_continuous(low = "beige", high = "firebrick") +
#    facet_wrap(~module_ix, scale = "free_y", ncol = 1) +
#    theme_minimal() +
#    theme(axis.text.y.right = element_text(size = 7, face = "italic", 
#					   margin = margin(l = .5, unit = "lines")),
#	  panel.grid = element_blank(),
#	  strip.text = element_blank(),
#	  plot.margin = margin(0, 0, 0, 0),
#	  strip.clip = "off"
#	  ) +
#    labs(x = NULL, y = NULL) +
#    coord_cartesian(clip = "off")


# Fig E #######################################################################

go_res <-
    "../bcell_lowinput/wgcna/data/DN2_go.tsv" |>
    read_tsv() |>
    group_by(module) |>
    top_n(5, -log10(pvalue)) |>
    ungroup() |>
    mutate(Description = str_trunc(Description, width = 30),
	   gene_r = map_dbl(GeneRatio, ~eval(parse(text = .)))) |>
    left_join(module_sizes, join_by(module)) |>
    select(module = ix, Description, pvalue, gene_r) |>
    arrange(module, pvalue) |>
    complete(module, fill = list(Description = NA)) |>
    mutate(Description = ifelse(is.na(Description), "No enrichment", Description),
	   Description = fct_inorder(Description))

fig_e <- 
    ggplot(data = go_res, 
       aes(x = " ", y = Description)) +
    geom_point(aes(fill = -log10(pvalue), size = gene_r), 
	       shape = 21, stroke = .5) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(position = "right") +
    scale_fill_gradient(low = "white", high = "black") +
    scale_size(range = c(.5, 3)) +
    ggforce::facet_col(vars(module), scales = "free_y", space = "free") + 
    theme_minimal() +
    theme(
	  axis.text.y.right = element_text(size = 8, margin = margin(0, .22, 0, 0, unit = "lines")),
	  axis.title = element_blank(),
	  legend.text = element_text(size = 8),
	  legend.title = element_text(size = 8),
	  legend.key.size = unit(.25, "cm"),
	  legend.margin = margin(r = 0.1, l = 0.1, unit = "lines"),
	  panel.grid = element_blank(),
	  plot.margin = margin(.5, 0, 0, 0, unit = "lines"),
	  strip.clip = "off",
	  strip.text = element_text(size = 8.5, margin = margin(b = 0))
	  ) +
    guides(fill = guide_colorbar("logP:", barheight = 5, barwidth = .25),
	   size = guide_legend("Gene\nRatio:"))


# Fig F #######################################################################
kme_df <- 
    "../bcell_lowinput/wgcna/data/DN2_kme.tsv" |>
    read_tsv() |>
    select(-grey) |>
    left_join(distinct(edger_res, gene_id, gene_name), join_by(gene_id)) |>
    pivot_longer(-c(gene_id, gene_name), names_to = "module", values_to = "kme") 

aid_genes <-
    tribble(~module, ~gene_name,
	    "black", "TNFSF4",
	    "black", "CDK2",
	    "blue", "TRAF1",
	    "blue", "IL2RA",
	    "brown", "BLK",
	    "brown", "IKBKE",
	    "green", "BCL11A",
	    "green", "BANK1",
	    "red", "SH2B3",
	    "red", "IL6R",
	    "yellow", "UBE2L3",
	    "yellow", "STAT1",
	    "turquoise", "ETS1",
	    "turquoise", "IL4R") |>
    left_join(module_sizes, join_by(module)) |>
    select(module, ix, gene_name) |>
    arrange(ix, gene_name)

aid_df <- 
    cpm_df |>
    filter(stim == "DN2") |>
    inner_join(aid_genes, join_by(gene_name)) |>
    select(sample_id, stim, condition, hours, gene_name, module = ix, cpm = obs_logcpm)
    
fig_f <-
    aid_df |>
    {function(x) split(x, x$module)}() |>
    {function(x) x[sort(names(x))]}() |>
    map(function(x) {

	    ggplot(data = x, 
		   aes(x = hours, y = cpm, fill = condition)) +
	    geom_quasirandom(aes(fill = condition),
			     method = "smiley", width = .2, 
			     shape = 21, stroke = .2, size = 1.5, alpha = .4) +
	    geom_line(aes(group = stim, color = stim), 
		      stat = "smooth", method = "loess", span = 1, se = FALSE, 
		      linewidth = .7) +
	    scale_y_continuous(breaks = scales::breaks_pretty(2)) +
	    scale_color_manual(values = pathway_colors) + 
	    scale_fill_manual(values = stim_colors) +
	    facet_wrap(~gene_name, nrow = 1, scales = "free_y") +
	    theme_minimal() +
	    theme(axis.text.x = element_blank(),
		  axis.text.y = element_text(size = 7, color = "grey70"),
		  panel.grid.major.y = element_blank(),
		  panel.grid.minor.y = element_blank(),
		  strip.text.x = element_text(size = 9, margin = margin(0, 0, .1, 0, unit = "lines")),
		  plot.title = element_text(size = 9, hjust = .5, margin = margin(0, 0, 0, 0)),
		  plot.margin = margin(.2, 0, 0, 0, unit = "lines"),
		  ) +
	    coord_cartesian(clip = "off") +
	    guides(fill = "none", color = "none") +
	    labs(x = NULL, y = NULL, title = unique(x$module))
	}
    ) |>
    {function(x) plot_grid(plotlist = x, ncol = 1)}() +
    theme(plot.margin = margin(t = .5, unit = "lines"))

bottom_title <- 
    ggdraw() + 
    draw_label(
	       "Gene programs in the DN2 condition: (d) expression patterns, (e) GO biological processes, and\n(f) example disease genes",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(t = 0, b = 0, l = 1.25, unit = "lines"))

bottom_fig <- 
    plot_grid(fig_d, NULL, fig_e, NULL, fig_f,
	      labels = c("d", "", "e", "", "f"),
	      label_size = 12, label_y = 1.05,
	      rel_widths = c(.66, .1, 1, .1, .66), nrow = 1) +
    theme(plot.margin = margin(t = .5, unit = "lines"))

bottom_grid <-
    plot_grid(bottom_title, bottom_fig, ncol = 1, rel_heights = c(.1, 1))


# Final fig
ggsave("./fig2.png",
       plot_grid(top_grid, NULL, fig_c, NULL, bottom_grid, 
		 ncol = 1, rel_heights = c(0.75, 0.05, 0.3, 0.025, 1)) +
       theme(plot.background = element_rect(fill = "white", color = "white")),
       width = 6.5, height = 9.2, dpi = 600)




#ldsc_df <- 
#    read_tsv("../bcell_lowinput/s-ldsc/compiled_results.tsv") |>
#    left_join(filter(module_sizes, stim == "DN2") |> select(module, ix), join_by(module))
#
#fig_f <- 
#    ggplot(data = ldsc_df,
#	   aes(x = factor(ix), y = gwas)) +
#    geom_tile(aes(fill = tau_star)) +
#    geom_point(data = filter(ldsc_df, pfdr <= 0.1),
#	       aes(x = ix, y = gwas),
#	       size = 1, shape = 8, color = "black") +
#    scale_y_discrete(labels = scales::label_wrap(20)) +
#    scale_fill_gradient2(low = "Navy Blue", mid = "white", high = "firebrick") +
#    theme_minimal() +
#    theme(
#	  axis.text = element_text(size = 8),
#	  axis.title = element_text(size = 8),
#	  axis.title.y = element_blank(),
#	  panel.grid = element_blank(),
#	  legend.text = element_text(size = 8),
#	  legend.title = element_text(size = 8),
#	  legend.position = "top",
#	  legend.margin = margin(l = -2.5, unit = "lines"),
#	  plot.margin = margin(1, 0, 0, 0, unit = "lines")) +
#    labs(x = "Module") +
#    guides(fill = guide_colorbar("Tau*:", barheight = .5, barwidth = 7, 
#				 title.position = "top"))
#
