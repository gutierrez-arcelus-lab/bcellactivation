# ==============================================================================
# Description: Generates Figure 2 (Low-input RNA-seq dynamics).
#              Panel A: PCA of samples over time.
#              Panel B: Heatmap of differentially expressed genes.
#              Panel C: B cell activation gene timecourses (cpm).
#              Panel D: WGCNA eigengene expression trajectories.
#              Panel E: GO biological processes enriched in modules.
#              Panel F: VST-normalized counts of WGCNA hub genes.
# ==============================================================================

library(tidyverse)
library(tidytext)
library(DESeq2)
library(cowplot)
library(ggbeeswarm)
library(glue)
library(ggh4x)
library(conflicted)
library(ggrastr)

conflicts_prefer(dplyr::filter())
conflicts_prefer(dplyr::slice())
conflicts_prefer(dplyr::count())
conflicts_prefer(dplyr::select())

# -----------------------------------------------------------------------------
# GLOBAL AESTHETICS & SETTINGS
# -----------------------------------------------------------------------------

# Set global ggplot theme to enforce max 7pt font size
theme_set(theme_minimal(base_size = 7))

stim_colors <-
    "./figure_colors_final.txt" |>
    read_tsv() |>
    mutate(lab = paste0(lab, "hrs")) |>
    select(lab, Hex) |>
    deframe()

stim_colors_0_names <- 
    keep(names(stim_colors), !grepl("Unstim", names(stim_colors))) |>
    str_remove("_\\d+hrs$") |>
    unique() |>
    paste0("_0hrs")

stim_colors_0 <- rep(stim_colors[[1]], length(stim_colors_0_names))
names(stim_colors_0) <- stim_colors_0_names

stim_colors <- c(stim_colors, stim_colors_0)
    
# -----------------------------------------------------------------------------
# Fig A: PCA
# -----------------------------------------------------------------------------
txi <- 
    "../01_rnaseq_lowinput/1_quantification/results/gene_expression_txi.rds" |>
    read_rds()

sample_table <- 
    tibble(sample_id = colnames(txi$counts)) |>
    separate(sample_id, c("donor_id", "stim", "time"), sep = "_", remove = FALSE) |>
    unite("condition", c(stim, time), sep = "_") |>
    column_to_rownames("sample_id")

dds <- 
    DESeqDataSetFromTximport(txi, sample_table, ~1) |>
    estimateSizeFactors()

pca_res <-
    plotPCA(vst(dds), ntop = 2000, returnData = TRUE)

pca_df <- 
    as_tibble(pca_res) |>
    select(sample_id = name, condition, PC1, PC2) |>
    separate(condition, c("stim", "time"), sep = "_", remove = FALSE) |>
    mutate(stim = recode(stim, 
			 "IL4" = "IL-4c", 
			 "CD40L" = "CD40c", 
			 "TLR9" = "TLR9c", 
			 "TLR7" = "TLR7c",
			 "BCR" = "BCRc", 
			 "BCR-TLR7" = "BCR/TLR7c", 
			 "DN2" = "DN2c")) |>
    unite("condition", c(stim, time), sep = "_", remove = FALSE) |>
    mutate(time = factor(time, levels = paste0(c(0, 4, 24, 48, 72), "hrs")),
	   condition = factor(condition, levels = names(stim_colors)))

pca_var <- 
    attr(pca_res, "percentVar") |>
    round(2) |>
    scales::percent()

legend_df <- 
    stim_colors |>
    enframe("condition", "color") |>
    separate("condition", c("stim", "time"), sep = "_", remove = FALSE) |>
    filter(condition %in% pca_df$condition) |>
    mutate(time = str_remove(time, "hrs"),
	   time = fct_inorder(time),
	   stim = fct_inorder(stim),
	   stim = fct_rev(stim)) |>
    select(stim, time, condition) |>
    mutate(stim = recode(stim, 
			 "IL4" = "IL-4c", 
			 "CD40L" = "CD40c", 
			 "TLR9" = "TLR9c", 
			 "TLR7" = "TLR7c",
			 "BCR" = "BCRc",
			 "BCR-TLR7" = "BCR/TLR7c",
			 "DN2" = "DN2c"),
	   stim = paste0(stim, ":"),
	   stim = fct_inorder(stim))

legend_plot <-
    ggplot(legend_df, aes(x = time, y = fct_rev(stim))) +
    geom_point(aes(fill = condition, shape = time), 
	       size = 2, stroke = .25) +
    scale_fill_manual(values = stim_colors) +
    scale_shape_manual(values = c(23, 25, 22, 21, 24)) +
    theme(
          axis.text.x = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.text.y = element_text(size = 7, margin = margin(r = -.5, unit = "lines")),
          axis.title.y =  element_blank(),
          axis.ticks.x = element_line(linewidth = .25),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(color = NA, fill = NA),
          plot.margin = margin(0, 0, 0, 0, unit = "pt"),
          plot.title = element_text(size = 7, margin = margin(t = 0, b = 0))
	) +
    labs(x = "Hours", title = "Conditions:") +
    guides(fill = "none", shape = "none")

pca_plot <- 
    ggplot(pca_df, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, color = "grey94") +
    geom_vline(xintercept = 0, color = "grey94") +
    stat_ellipse(aes(group = time), level = .85,
		 linetype = 2, linewidth = 0.25, color = "grey50") +
    geom_text(data = 
	       tribble(~x, ~y, ~lab,
		       0, 50, "0h",
		       -75, 0, "4h",
		       -25, -50, "24h",
		       25, -75, "48h",
		       42, -50, "72h"),
	       aes(x, y, label = lab),
	       size = 7, size.unit = "pt") +
    geom_point(aes(fill = condition, shape = time), 
	       size = 2, stroke = .2) +
    scale_fill_manual(values = stim_colors) +
    scale_shape_manual(values = c(23, 25, 22, 21, 24)) +
    theme(
	  axis.text = element_blank(),
	  axis.title = element_text(size = 7),
	  panel.grid = element_blank(),
	  plot.margin = margin(0, 0, 0, 0, unit = "lines")) +
    guides(fill = "none", shape = "none") +
    labs(x = glue("PC1 ({pca_var[1]})"), y = glue("PC2 ({pca_var[2]})"))

pca_inset <- 
    ggplot(pca_df |> mutate(i = stim %in% c("BCRc", "BCR/TLR7c", "DN2c")),
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
	  plot.margin = margin(t = 0, r = 1.25, b = 0, l = 1, unit = "lines")) +
    guides(color = "none") +
    labs(title = "BCR stim:")

fig_a_title <- 
    ggdraw() + 
    draw_label(
	       "PCA shows separation of stimuli and time points\n ",
	       x = 0,
	       size = 7,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 7),
	  plot.margin = margin(t = 0.5, l = 1.25, unit = "lines"))


# -----------------------------------------------------------------------------
# Fig B: Differential Expression Heatmap
# -----------------------------------------------------------------------------
all_stims <-
    names(stim_colors) |>
    str_split("_") |>
    map_chr(1) |>
    unique()

diff_expr <- 
    "../01_rnaseq_lowinput/2_dge/results/edger/diff_expr_all_times.tsv" |>
    read_tsv()

diff_expr_summ <-
    diff_expr |>
    group_by(group1, group2) |>
    summarise(n = sum(!is.na(gene_id))) |>
    ungroup() |>
    separate(group1, c("stim1", "t1"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    separate(group2, c("stim2", "t2"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    mutate_at(vars(stim1, stim2), ~recode(., "BCR_TLR7" = "BCR-TLR7")) |>
    mutate_at(vars(stim1, stim2), ~recode(., 
					 "IL4" = "IL-4c", 
					 "CD40L" = "CD40c", 
					 "TLR9" = "TLR9c", 
					 "TLR7" = "TLR7c",
					 "BCR" = "BCRc", 
					 "BCR-TLR7" = "BCR/TLR7c", 
					 "DN2" = "DN2c")) |>
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
    scico::scale_fill_scico(palette = "nuuk", 
			    na.value = "white",
			    labels = scales::comma) +
    theme(
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  panel.grid.major = element_blank(),
	  panel.background = element_rect(color = NA, fill = "transparent"),
	  plot.margin = margin(0, 0, 0, 0),
	  legend.text = element_text(size = 7),
	  legend.title = element_text(size = 7),
	  legend.position.inside = c(.8, .3),
	  ) +
    labs(fill = "DE\ngenes:") +
    guides(fill = guide_colorbar(barwidth = .5, barheight = 5, position = "inside")) +
    coord_cartesian(xlim = c(1, 28), ylim = c(1, 28))

b_axis_x_df <- 
    diff_expr_summ |>
    filter(!is.na(stim1)) |>
    distinct(group1, stim1, t1) |>
    mutate(condition1 = glue("{stim1}_{t1}hrs")) |>
    select(group1, condition1)

b_axis_y_df <- 
    diff_expr_summ |>
    filter(!is.na(stim2)) |>
    distinct(group2, stim2, t2) |>
    mutate(condition2 = glue("{stim2}_{t2}hrs")) |>
    select(group2, condition2)

b_axis_x_plot <-
    ggplot(data = b_axis_x_df, 
	   aes(y = factor(1), x = group1)) +
    geom_tile(aes(fill = condition1), color = "black") +
    scale_fill_manual(values = stim_colors) +
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
    geom_tile(aes(fill = condition2), color = "black") +
    scale_fill_manual(values = stim_colors) +
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
	       "Number of differentially expressed genes across all\nconditions and time points",
	       x = 0,
	       size = 7,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 7),
	  plot.margin = margin(l = 1.25, unit = "lines"))
		 
b_tmp <- 
    plot_grid(NULL, b_axis_x_plot, b_axis_y_plot, diff_plot,
	      ncol = 2, nrow = 2, align = "v", rel_heights = c(.075, 1), 
	      rel_widths = c(0.06, 1))


# -----------------------------------------------------------------------------
# Fig C: Gene Expression Timecourses
# -----------------------------------------------------------------------------
gene_exp <-
    "../01_rnaseq_lowinput/2_dge/results/edger/cpm.tsv" |>
    read_tsv() |>
    extract(sample_id, "timep", ".+_.+_(\\d+)hrs", convert = TRUE, remove = FALSE) |>
    select(stim, timep, sample_name = sample_id, gene_id, gene_name, cpm = obs_cpm)    

count_data <- 
    gene_exp |> 
    mutate(stim = recode(stim, 
                         "IL4" = "IL-4c",
                         "CD40L" = "CD40c",
                         "TLR9" = "TLR9c",
                         "TLR7" = "TLR7c",
                         "BCR" = "BCRc",
                         "BCR-TLR7" = "BCR/TLR7c",
                         "DN2" = "DN2c"),
           stim = factor(stim, 
                         levels = c("IL-4c", "CD40c", "TLR9c", "TLR7c",
                                    "BCRc", "BCR/TLR7c", "DN2c")),
           timep = factor(timep, levels = c("0", "4", "24", "48", "72")),
           condition = glue("{stim}_{timep}hrs")
    ) |>
    arrange(stim, timep)

pathway_colors <-
    keep(stim_colors, grepl("48hrs", names(stim_colors)))

names(pathway_colors) <- str_remove(names(pathway_colors), "_48hrs")

pathway_colors <- c("IL-4c" = "goldenrod4", pathway_colors)

timecourse_plot_1 <- 
    ggplot(data = count_data |> filter(gene_name %in% c("AICDA", "FCER2", "BANK1")), 
           aes(x = timep, y = cpm)) +
    geom_quasirandom(aes(fill = condition),
                     method = "smiley", width = .2, 
                     shape = 21, stroke = .1, size = 1.25, alpha = .4) +
    geom_line(aes(group = stim, color = stim), 
              stat = "smooth", method = "loess", span = 1, se = FALSE, 
              linewidth = .7) +
    scale_y_continuous(limits = ~ c(floor(min(.x)), ceiling(max(.x))),
                       breaks = ~ c(max(c(0, ceiling(.x[1]))), floor(.x[2]))) +
    scale_color_manual(values = pathway_colors) + 
    scale_fill_manual(values = stim_colors) + 
    facet_grid(gene_name~stim, scale = "free", space = "free_x") +
    theme(
	  axis.text.x = element_blank(),
	  axis.text.y = element_text(size = 7),
	  axis.title.y = element_text(size = 7),
	  strip.text.x = element_blank(),
	  strip.text.y = element_text(size = 7, angle = 0, face = "italic"),
	  strip.clip = "off",
	  panel.grid.minor = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.spacing.x = unit(0, "null"),
	  panel.spacing.y = unit(.5, "lines"),
	  legend.position = "none",
	  plot.margin = margin(t = 0, r = .5, b = 0, l = 0, unit = "lines")
	  ) + 
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Counts per Million")

timecourse_plot_2 <- 
    ggplot(data = count_data |> filter(gene_name %in% c("IRF5", "TBX21", "XBP1")), 
	   aes(x = timep, y = cpm)) +
    geom_quasirandom(aes(fill = condition),
                     method = "smiley", width = .2, 
                     shape = 21, stroke = .1, size = 1.25, alpha = .4) +
    geom_line(aes(group = stim, color = stim), 
              stat = "smooth", method = "loess", span = 1, se = FALSE, 
              linewidth = .7) +
    scale_y_continuous(limits = ~ c(floor(min(.x)), ceiling(max(.x))),
                       breaks = ~ c(max(c(0, ceiling(.x[1]))), floor(.x[2]))) +
    scale_color_manual(values = pathway_colors) + 
    scale_fill_manual(values = stim_colors) + 
    facet_grid(gene_name~stim, scale = "free", space = "free_x") +
    theme(
	  axis.text.x = element_blank(),
	  axis.text.y = element_text(size = 7),
	  axis.title.y = element_text(size = 7),
	  strip.text.x = element_blank(),
	  strip.text.y = element_text(size = 7, angle = 0, face = "italic"),
	  strip.clip = "off",
	  panel.grid.minor = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.spacing.x = unit(0, "null"),
	  panel.spacing.y = unit(.5, "lines"),
	  legend.position = "none", 
	  plot.margin = margin(t = 0, r = 0, b = 0, l = 0.5, unit = "lines")
    ) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Counts per Million")

fig_c_grid <- 
    plot_grid(timecourse_plot_1 , timecourse_plot_2, nrow = 1)

fig_c_title <- 
    ggdraw() + 
    draw_label(
	       "B cell activation genes show distinct patterns of expression across conditions",
	       x = 0,
	       size = 7,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 7),
	  plot.margin = margin(b = .05, l = 1.25, unit = "lines"))
       

# -----------------------------------------------------------------------------
# Fig D: WGCNA Modules
# -----------------------------------------------------------------------------
module_sizes <-
    read_tsv("../01_rnaseq_lowinput/3_wgcna/data/DN2_modules.tsv") |>
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
    "../01_rnaseq_lowinput/3_wgcna/data/DN2_eigen.tsv" |>
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

kim_df <- 
    "../01_rnaseq_lowinput/3_wgcna/data/DN2_kim.tsv" |>
    read_tsv() |>
    left_join(module_sizes, join_by(module)) |> 
    select(gene_id, gene_name, module, module_ix = ix, kim) |>
    arrange(module_ix, desc(kim)) |>
    group_by(module) |>
    top_n(5, kim) |>
    ungroup()

module_plot <-
    ggplot(modules_df) +
    geom_line(aes(x = time, y = value, group = module),
	      linewidth = 1.25) +
    facet_wrap(~module_ix, ncol = 1) +
    theme(
	  axis.text.x = element_text(size = 7),
	  axis.text.y = element_text(size = 7),
	  axis.title.x = element_blank(),
	  axis.title.y = element_text(size = 7),
	  panel.grid.major.x = element_line(linetype = 2, color = "black", linewidth = .1),
	  panel.grid.major.y = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  panel.spacing.y = unit(0.15, "lines"),
	  strip.text = element_text(size = 7, margin = margin(t = 0, b = 0)),
	  strip.clip = "off",
	  plot.margin = margin(0, 0, 0, 0, unit = "lines")) +
    labs(y = "Eigengene expression")

kim_plot <-
    ggplot(data = kim_df |> 
	   arrange(module_ix, kim) |>
	   mutate(gene_name = fct_inorder(gene_name)),
	   aes(x = " ", y = gene_name)) +
    geom_point(aes(fill = kim), shape = 21, size = 7/.pt, show.legend = FALSE) +
    scale_x_discrete(drop = FALSE, expand = c(0, 0), limits = c(" ", " ")) +
    scale_y_discrete(position = "right", expand = c(0, 0)) +
    scale_fill_continuous(low = "beige", high = "firebrick") +
    facet_wrap(~module_ix, scale = "free_y", ncol = 1,
	       labeller = as_labeller(c("Module 1" = " ", "Module 2" = " ", "Module 3" = " ",
					"Module 4" = " ", "Module 5" = " ", "Module 6" = " ",
					"Module 7" = " ", "Module 8" = " "))) +
    theme(
	  axis.text.y.right = element_text(size = 7, face = "italic", 
					   margin = margin(l = .5, unit = "lines")),
	  panel.grid = element_blank(),
	  panel.spacing.y = unit(0.15, "lines"),
	  plot.margin = margin(0, 0, 0, 0, unit = "lines"),
	  strip.text = element_text(size = 8.5, margin = margin(t = 0, b = 0)),
	  strip.clip = "off"
	  ) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off")

fig_d <- 
    plot_grid(module_plot, NULL, kim_plot, nrow = 1, rel_widths = c(1, .1, .5))

# -----------------------------------------------------------------------------
# Fig E: GO Enrichment
# -----------------------------------------------------------------------------
go_res <-
    "../01_rnaseq_lowinput/3_wgcna/data/DN2_go.tsv" |>
    read_tsv() |>
    group_by(module) |>
    top_n(5, -log10(pvalue)) |>
    ungroup() |>
    mutate(Description = str_trunc(Description, width = 36),
	   gene_r = map_dbl(GeneRatio, ~eval(parse(text = .)))) |>
    right_join(module_sizes, join_by(module)) |>
    select(module = ix, Description, pvalue, gene_r) |>
    arrange(module, pvalue) |>
    mutate(Description = fct_inorder(Description))

fig_e <- 
    ggplot(data = go_res, 
       aes(x = " ", y = Description)) +
    geom_point(aes(fill = -log10(pvalue), size = gene_r), 
	       shape = 21, stroke = .5) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(position = "right") +
    scale_fill_gradient(low = "white", high = "black") +
    scale_size(range = c(.5, 2.5)) +
    ggforce::facet_col(vars(module), scales = "free_y", drop = FALSE) + 
    theme(
	  text = element_text(size = 7),
	  axis.text.y.right = element_text(size = 7, margin = margin(r = .5, l = 0, unit = "lines")),
	  axis.title = element_blank(),
	  legend.text = element_text(size = 7),
	  legend.title = element_text(size = 7),
	  legend.key.size = unit(.25, "cm"),
	  legend.margin = margin(r = 0.1, l = 0.1, unit = "lines"),
	  legend.box.spacing = unit(1, "pt"),
	  panel.grid = element_blank(),
	  plot.margin = margin(0, 0, 0, 0, unit = "lines"),
	  strip.clip = "off",
	  strip.text = element_text(size = 7, margin = margin(t = 0, b = 0, l = 2, unit = "lines")),
	  panel.spacing.y = unit(0.15, "lines")
	  ) +
    guides(fill = guide_colorbar("logP:", barheight = 5, barwidth = .25),
	   size = guide_legend("Gene\nRatio:"))

# -----------------------------------------------------------------------------
# Fig F: WGCNA Hub Gene Examples (RESTORED TO NESTED FACETS)
# -----------------------------------------------------------------------------
gene_metadata <-
    count_data |>
    distinct(gene_id, gene_name)

kme_df <- 
    "../01_rnaseq_lowinput/3_wgcna/data/DN2_kme.tsv" |>
    read_tsv() |>
    select(-grey) |>
    left_join(gene_metadata, join_by(gene_id)) |>
    pivot_longer(-c(gene_id, gene_name), names_to = "module", values_to = "kme") |>
    group_by(gene_id) |>
    slice_max(kme) |>
    ungroup()

aid_genes <-
    tribble(~gene_name,
        "TNFSF4", "IL12RB1", "TRAF1", "IL2RA", "CSK", "IKBKE",
        "BLK", "BANK1", "SH2B3", "NFKBIA", "UBE2L3", "IRF5",
        "ETS1", "IL4R", "ITGA4", "CXCR4") |>
    inner_join(kme_df, join_by(gene_name)) |>
    left_join(module_sizes, join_by(module)) |>
    select(module, ix, gene_name) |>
    arrange(ix, gene_name)

vst_df <- 
    vst(dds) |>
    assay() |> 
    as.data.frame() |>
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "sample_name", values_to = "vst_counts") |>
    select(sample_name, gene_id, vst_counts)

aid_df <- 
    count_data |>
    select(-cpm) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    inner_join(vst_df, join_by(sample_name, gene_id)) |>
    inner_join(aid_genes, join_by(gene_name)) |>
    filter(stim == "DN2c") |>
    select(sample_name, condition, stim, timep, gene_id, gene_name, module = ix, vst_counts)

# 1. Pre-compute the 16 unique y-axis scales sequentially
y_scales <- 
    aid_df |>
    arrange(module, gene_name) |>
    group_by(module, gene_name) |>
    summarise(mn = min(vst_counts), mx = max(vst_counts), .groups = "drop") |>
    mutate(scale_obj = map2(mn, mx, ~ scale_y_continuous(
                       expand = expansion(mult = c(0.1, 0.1)),
                       limits = c(floor(.x), max(ceiling(.y), 0.1)),
                       breaks = c(max(floor(.x), 0), ceiling(.y))))) |>
    pull(scale_obj)

# 2. Build the nested strip configuration
strips <- strip_nested(
    text_x = list(element_text(size = 7), element_text(size = 7, face = "italic")),
    by_layer_x = TRUE
)

fig_f <- 
    ggplot(data = aid_df, aes(x = timep, y = vst_counts)) +
    geom_quasirandom(aes(fill = condition),
		     method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 1.5, alpha = .4) +
    geom_line(aes(group = stim, color = stim), 
	      stat = "smooth", method = "loess", span = 1, se = FALSE, 
	      linewidth = .7) +
    scale_color_manual(values = pathway_colors) + 
    scale_fill_manual(values = stim_colors) +
    facet_nested_wrap(~ module + gene_name, 
              scales = "free_y", ncol = 2,
              strip = strips) +
    facetted_pos_scales(y = y_scales) +
    theme(
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7),
          axis.title = element_text(size = 7),
	  panel.grid.major.x = element_line(linetype = 2, color = "black", linewidth = .1),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
	  panel.spacing.y = unit(0.15, "lines"),
          strip.text = element_text(size = 7, margin = margin(t = 0, b = 0)),
          plot.margin = margin(t = 0, 0, 0, 0, unit = "lines")
          ) +
    coord_cartesian(clip = "off") +
    guides(fill = "none", color = "none") +
    labs(x = "Time", y = "VST-normalized counts")

# -----------------------------------------------------------------------------
# Compilation & Export
# -----------------------------------------------------------------------------
bottom_title <- 
    ggdraw() + 
    draw_label(
	       "Gene programs in DN2c: (d) expression patterns and hub genes, (e) GO biological processes, and (f) example disease genes",
	       x = 0,
	       size = 7,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 7),
	  plot.margin = margin(t = 0, b = .5, l = 1.25, unit = "lines"))

# 1. Apply the 1.25 line left-margin directly to the internal plots/grids
plots_a_indented <- 
    plot_grid(pca_plot, 
	      plot_grid(pca_inset, NULL, legend_plot, ncol = 1, rel_heights = c(.4, .05, 1)), 
	      nrow = 1, rel_widths = c(1, .5)) + 
    theme(plot.margin = margin(l = 0.5, unit = "lines"))

b_tmp_indented <- 
    b_tmp + 
    theme(plot.margin = margin(t = 0.5, r = 0.5, unit = "lines"))

fig_c_indented <- 
    fig_c_grid + 
    theme(plot.margin = margin(t = 0.5, r = 0.5, b = .5, l = 0.5, unit = "lines"))

fig_d_indented <- 
    fig_d + 
    theme(plot.margin = margin(0, 0, 0, 0.5, unit = "lines"))

fig_e_indented <- 
    fig_e + 
    theme(plot.margin = margin(0, 0, 0, 0, unit = "lines"))

fig_f_indented <- 
    fig_f + 
    theme(plot.margin = margin(0, 0.5, 0, 0, unit = "lines"))

# 2. Rebuild the labeled grids over the indented plots (Labels sit at x=0)
fig_a_grid <- 
    plot_grid(fig_a_title, plots_a_indented, 
	      ncol = 1, rel_heights = c(.1, 1), 
	      labels = c("a", NULL), label_size = 10, label_y = .9, vjust = 1)

fig_b_grid <- 
    plot_grid(fig_b_title, b_tmp_indented, 
	      ncol = 1, rel_heights = c(.1, 1), 
	      labels = c("b"), label_size = 10, label_y = .9, vjust = 1)

fig_c <- 
    plot_grid(fig_c_title, fig_c_indented, 
	      ncol = 1, rel_heights = c(.1, 1), 
	      labels = "c", label_size = 10, label_y = .9, vjust = 1)

# 3. Assemble the main rows
top_grid <- 
    plot_grid(fig_a_grid, NULL, fig_b_grid, 
	      nrow = 1, rel_widths = c(1, .05, .8))

bottom_fig <-
    plot_grid(fig_d_indented, NULL, fig_e_indented, NULL, fig_f_indented,
	      labels = c("d", "", "e", "", "f"), label_size = 10, label_y = 1, vjust = 1, 
	      rel_widths = c(.8, .025, 1, .1, .7), nrow = 1)

bottom_grid <- 
    plot_grid(bottom_title, bottom_fig, ncol = 1, rel_heights = c(.025, 1))

# 4. Final Combination
final_figure <-
    plot_grid(top_grid, NULL, fig_c, NULL, bottom_grid, 
	      ncol = 1, rel_heights = c(0.55, 0.03, 0.25, 0.01, 1))

# Final fig
ggsave("./pdf/fig2.pdf", 
       final_figure,
       width = 179,
       height = 220,
       units = "mm",
       dpi = 600,
       device = cairo_pdf
       )
