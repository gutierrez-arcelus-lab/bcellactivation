library(tidyverse)
library(tidytext)
library(DESeq2)
library(extrafont)
library(furrr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(cowplot)
library(ggbeeswarm)
library(glue)
library(ggraph)
library(tidygraph)

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
    select(stim, time, condition) |>
    mutate(stim = paste0(stim, ":"),
	   stim = fct_inorder(stim))

legend_plot <-
    ggplot(legend_df, aes(x = time, y = 1, fill = condition)) +
    geom_point(size = 2.5, shape = 21, stroke = .25) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 1) +
    theme_minimal() +
    theme(text = element_text(size = 9),
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
	  plot.margin = margin(l = 0, r = 0, unit = "cm")
	  ) +
    labs(x = "hours") +
    guides(fill = "none")

pca_plot <- 
    ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(fill = condition), 
	       size = 2.5, shape = 21, stroke = .15) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  panel.grid = element_blank(),
	  plot.margin = margin(r = 0),
	  plot.title = element_text(size = 9)) +
    guides(fill = "none") +
    labs(x = pc_varexp$lab[1], y = pc_varexp$lab[2])

fig_a_title <- 
    ggdraw() + 
    draw_label(
	       "PCA shows separation of stimuli and time points",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1, unit = "lines"))

fig_a_grid <- 
    plot_grid(
	      fig_a_title,
	      plot_grid(
			pca_plot, 
			plot_grid(NULL, legend_plot, NULL, ncol = 1, rel_heights = c(.1, 1, .1)), 
			nrow = 1, 
			rel_widths = c(1, .33)
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

diff_plot <- 
    ggplot(data = diff_expr_summ, 
       aes(x = group1, y = group2)) +
    geom_tile(aes(fill = n), alpha = 1) +
    scale_fill_viridis_c(option = "cividis",
			 na.value = "white",
			 labels = scales::comma) +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  panel.grid.major = element_blank(),
	  panel.background = element_rect(fill = "transparent", color = NA),
	  legend.text = element_text(margin = margin(l = 0)),
	  legend.box.spacing = unit(c(l = -.25), "lines"),
	  legend.margin = margin(0, 0, 0, -.25, unit = "lines"),
	  plot.margin = margin(0, 0, 0, 0)
	  ) +
    labs(fill = "DE\ngenes:") +
    guides(fill = guide_colorbar(barwidth = .5, barheight = 7)) +
    coord_cartesian(xlim = c(1, 28), ylim = c(1, 28))

b_legend <- get_plot_component(diff_plot, 'guide-box-right', return_all = TRUE)

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
	  plot.margin = margin(b = .5, l = 1.25, unit = "lines"))
		 
b_tmp <- 
    plot_grid(NULL, b_axis_x_plot, b_axis_y_plot, diff_plot + guides(fill = "none"),
	      ncol = 2, nrow = 2, align = "v", rel_heights = c(.05, 1), rel_widths = c(0.05, 1))

fig_b_grid <- 
    plot_grid(
	      fig_b_title, 
	      plot_grid(b_tmp, b_legend, rel_widths = c(1, .175), nrow = 1),
	      ncol = 1, 
	      rel_heights = c(.15, 1),
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
    scale_y_continuous(expand = expansion(mult = c(0.25, 0)),
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
	      plot_grid(plot_timecourse(filter(cpm_df, gene_name == "BANK1")),
			plot_timecourse(filter(cpm_df, gene_name == "IRAK1")),
			plot_timecourse(filter(cpm_df, gene_name == "IKZF1")),
			ncol = 1),
	      NULL,
	      plot_grid(plot_timecourse(filter(cpm_df, gene_name == "STAT1")),
			plot_timecourse(filter(cpm_df, gene_name == "TBX21")),
			plot_timecourse(filter(cpm_df, gene_name == "TNFSF4")),
			ncol = 1),
	      NULL,
	      nrow = 1, rel_widths = c(0.05, 1, 0.05, 1, 0.05)) +
    draw_label("Counts per million", x = -.001, y = 0.5, vjust = 1.5, angle = 90, size = 9)
       
fig_c <- 
    plot_grid(fig_c_title, fig_c_grid, ncol = 1, rel_heights = c(.1, 1),
	      labels = "c", label_size = 12)




# Fig D #######################################################################
module_sizes <-
    list.files("../bcell_lowinput/wgcna/data",
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

module_colors <- unique(module_sizes$module)
module_colors <- setNames(module_colors, module_colors)
module_colors["yellow"] <- "goldenrod2"
module_colors["green"] <- "forestgreen"

eigengenes_df <-
    list.files("../bcell_lowinput/wgcna/data",
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

dn2_modules_df <-
    eigengenes_df |>
    filter(stim == "DN2") |>
    group_by(stim, module, module_ix, time) |>
    summarise(value = mean(value)) |>
    ungroup() |>
    arrange(stim, module_ix, time) |>
    mutate(module_ix = fct_inorder(as.character(module_ix)))

kim_df <- 
    "../bcell_lowinput/wgcna/data/DN2_kim.tsv" |>
    read_tsv() |>
    left_join(filter(module_sizes, stim == "DN2"), join_by(module)) |> 
    select(gene_id, gene_name, module, module_ix = ix, kim) |>
    arrange(module_ix, desc(kim)) |>
    group_by(module) |>
    top_n(5, kim) |>
    ungroup()

modules_plot <-
    ggplot(dn2_modules_df) +
    geom_line(aes(x = time, y = value, group = module, color = module),
	      linewidth = 1.25) +
    scale_color_manual(values = module_colors) +
    facet_wrap(~module_ix, ncol = 1, labeller = labeller(module_ix = function(i) paste("Module", i))) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8),
	  axis.text.y = element_text(size = 8, margin = margin(r = -1, l = 1, unit = "lines")),
	  axis.title.x = element_blank(),
	  axis.title.y = element_text(size = 8, margin = margin(r = -.5, l = .5, unit = "lines")),
	  panel.grid.major.x = element_line(linetype = 2, color = "black", linewidth = .1),
	  panel.grid.major.y = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  strip.text = element_text(size = 8.5, margin = margin(b = 0)),
	  strip.clip = "off",
	  plot.margin = margin(1, 0, 0, 0, unit = "lines")) +
    guides(color = "none") +
    labs(y = "Eigengene expression")


kim_plot <-
    ggplot(data = kim_df |> 
	   arrange(module_ix, kim) |>
	   mutate(gene_name = fct_inorder(gene_name)),
	   aes(x = " ", y = gene_name)) +
    geom_point(aes(fill = kim), shape = 21, size = 3, show.legend = FALSE) +
    scale_x_discrete(drop = FALSE, expand = c(0, 0), limits = c(" ", " ")) +
    scale_y_discrete(position = "right") +
    scale_fill_continuous(low = "beige", high = "firebrick") +
    facet_wrap(~module_ix, scale = "free_y", ncol = 1) +
    theme_minimal() +
    theme(axis.text.y.right = element_text(size = 7, face = "italic", 
					   margin = margin(l = .5, unit = "lines")),
	  panel.grid = element_blank(),
	  strip.text = element_blank(),
	  plot.margin = margin(1, 0, 0, 0, unit = "lines"),
	  strip.clip = "off"
	  ) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off")

fig_d <-  
    plot_grid(modules_plot, NULL, kim_plot, 
	      nrow = 1, rel_widths = c(1, .1, .6))


# Fig E #######################################################################
go_res <-
    "../bcell_lowinput/wgcna/data/DN2_go.tsv" |>
    read_tsv() |>
    group_by(module) |>
    top_n(6, -log10(pvalue)) |>
    ungroup() |>
    mutate(stim = "DN2",
	   Description = str_trunc(Description, width = 25),
	   gene_r = map_dbl(GeneRatio, ~eval(parse(text = .)))) |>
    left_join(module_sizes, join_by(stim, module)) |>
    select(module = ix, Description, pvalue, gene_r) |>
    mutate(module = factor(module, levels = sort(unique(module_sizes$ix)))) 

fig_e <- 
    ggplot(data = go_res, 
       aes(x = " ", 
	   y = reorder_within(Description, by = -log10(pvalue), within = module))) +
    geom_point(aes(fill = -log10(pvalue), size = gene_r), 
	       shape = 21, stroke = .5) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_reordered(position = "right") +
    scale_fill_gradient(low = "white", high = "black") +
    scale_size(range = c(.5, 3)) +
    facet_wrap(~module, ncol = 1, scale = "free_y", 
	       labeller = labeller(module = function(i) paste("Module", i))) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
	  axis.text.y.right = element_text(size = 8, margin = margin(0, 0, 0, -2)),
	  legend.text = element_text(size = 8, margin = margin(l = -0, unit = "lines")),
	  legend.title = element_text(size = 8),
	  legend.key.size = unit(.25, "cm"),
	  panel.grid = element_blank(),
	  plot.margin = margin(1, 0, 0, 0, "lines"),
	  strip.clip = "off",
	  strip.text.x = element_text(size = 8.5, hjust = 1, margin = margin(.25, -6, 0, 0, "lines")),
	  ) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colorbar("logP:", barheight = 5, barwidth = .25),
	   size = guide_legend("Gene\nRatio:"))


# Fig F #######################################################################
ldsc_files <- 
    list.files("../bcell_lowinput/s-ldsc/results", 
	       pattern = "cell_type_results.txt",
	       full.names = TRUE)

names(ldsc_files) <- ldsc_files |> 
    basename() |>
    str_remove("modules_500_seg_") |>
    str_remove(".cell_type_results.txt")

ldsc_df <- 
    map_dfr(ldsc_files, read_tsv, .id = "gwas") |>
    mutate(pfdr = p.adjust(Coefficient_P_value, method = "fdr"),
	   gwas = str_remove(gwas, "PASS_"),
	   gwas = str_remove(gwas, "UKB_460K\\.(body|disease)_"),
	   gwas = str_replace_all(gwas, "_", " "),
	   gwas = str_to_title(gwas),
	   gwas = recode(gwas, "Heightz" = "Height")) |>
    select(gwas, module = Name, cf = Coefficient, pfdr) |>
    left_join(filter(module_sizes, stim == "DN2") |> select(module, ix), join_by(module))

fig_f <- 
    ggplot(data = ldsc_df,
	   aes(x = factor(ix), y = gwas)) +
    geom_tile(aes(fill = cf)) +
    geom_point(data = filter(ldsc_df, pfdr <= 0.1),
	       aes(x = ix, y = gwas),
	       size = 1, shape = 8, color = "black") +
    scale_y_discrete(labels = scales::label_wrap(20)) +
    scale_fill_gradient2(low = "Navy Blue", mid = "white", high = "firebrick") +
    theme_minimal() +
    theme(
	  axis.text = element_text(size = 8),
	  axis.title = element_text(size = 8),
	  axis.title.y = element_blank(),
	  panel.grid = element_blank(),
	  legend.text = element_text(size = 8),
	  legend.title = element_text(size = 8),
	  legend.position = "top",
	  legend.margin = margin(l = -2.5, unit = "lines"),
	  plot.margin = margin(1, 0, 0, 0, unit = "lines")) +
    labs(x = "Module") +
    guides(fill = guide_colorbar("Coefficient:", barheight = .5, barwidth = 7.5, 
				 title.position = "top"))

bottom_title <- 
    ggdraw() + 
    draw_label(
	       "Gene programs in the DN2 condition: (d) expression pattern and hub genes, (e) biological\nprocesses, and (f) enrichment of disease heritability.",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(t = 0, b = 0, l = 1.25, unit = "lines"))

bottom_fig <- 
    plot_grid(fig_d, NULL, fig_e, NULL, fig_f,
	      labels = c("d", " ", "e", " ", "f"),
	      label_size = 12,
	      rel_widths = c(1, .1, 1, .1, 1), nrow = 1)

bottom_grid <-
    plot_grid(bottom_title, bottom_fig, ncol = 1, rel_heights = c(.1, 1))


# Final fig
ggsave("./testp.png",
       plot_grid(top_grid, NULL, fig_c, NULL, bottom_grid, 
		 ncol = 1, rel_heights = c(0.7, 0.05, 0.45, 0.05, 1)) + 
       theme(plot.background = element_rect(fill = "white", color = "white")),
       width = 6.5, height = 9.5)


