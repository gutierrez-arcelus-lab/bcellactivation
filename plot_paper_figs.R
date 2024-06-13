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

slice <- dplyr::slice
count <- dplyr::count
select <- dplyr::select

if (!file.exists("paper_plots")) dir.create("paper_plots")



# Left-hand side ##############################################################

# Fig A PCA

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
	       size = 2.5, shape = 21, stroke = .15) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  panel.grid = element_blank(),
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
	  plot.margin = margin(0, 1.25, 0, 2, unit = "lines"))

fig_a_grid <- 
    plot_grid(
	      fig_a_title,
	      plot_grid(
			pca_plot, 
			plot_grid(NULL, legend_plot, NULL, ncol = 1, rel_heights = c(.1, 1, .5)), 
			nrow = 1, 
			rel_widths = c(1, .6)
			),
	      ncol = 1, 
	      rel_heights = c(.1, 1),
	      labels = c(NULL, "a"), label_size = 10
    )


# Fig C DN2 modules
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
    mutate(module = factor(module, levels = sort(unique(module_sizes$ix)))) 

dn2_modules_df <-
    eigengenes_df |>
    filter(stim == "DN2") |>
    filter(module_ix %in% unique(go_res$module)) |>
    group_by(stim, module = module_ix, time) |>
    summarise(value = mean(value)) |>
    ungroup() |>
    mutate(module = fct_inorder(as.character(module)))

dn2_modules_plot <- 
    ggplot() +
    geom_line(data = dn2_modules_df,
	      aes(x = time, y = value, group = module),
	      linewidth = 1.25) +
    facet_wrap(~module, ncol = 1) +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  axis.text = element_blank(),
	  axis.ticks.y = element_blank(),
	  panel.grid = element_blank(),
	  strip.text = element_text(size = 9,
				    margin = margin(t = 0, b = 0)),
	  strip.background = element_rect(fill = "grey80", color = NA),
	  plot.margin = unit(c(5.5, 0, 5.5, 5.5), "pt")
	  ) +
    labs(x = NULL, y = "Average Eigengene expression")


kim_df <- 
    "./bcell_lowinput/wgcna/data/DN2_kim.tsv" |>
    read_tsv() |>
    group_by(module) |>
    top_n(5, kim) |>
    ungroup() |>
    mutate(stim = "DN2") |>
    left_join(module_sizes, join_by(stim, module)) |>
    select(module = ix, gene_name, kim) |>
    filter(module %in% unique(go_res$module)) |>
    arrange(module, desc(kim))

kim_plot <-
    ggplot(data = kim_df,
	   aes(x = 1, 
	       y = reorder_within(gene_name, by = kim, within = module))) +
    geom_tile(aes(fill = kim), show.legend = FALSE) +
    scale_y_reordered(position = "right") +
    scale_fill_continuous(low = "beige", high = "firebrick") +
    facet_wrap(~module, scale = "free_y", ncol = 1) +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  axis.text.x = element_blank(),
	  axis.text.y.right = element_text(hjust = 0, 
					   margin = margin(l = -.25, unit = "lines")),
	  panel.grid = element_blank(),
	  strip.text = element_blank(),
	  plot.margin = unit(c(5.5, 5.5, 5.5, 2.5), "pt")
	  ) +
    labs(x = NULL, y = NULL)

# GO
go_plot <- 
    ggplot(data = go_res, 
       aes(x = "1", 
	   y = reorder_within(Description, by = -log10(pvalue), within = module))) +
    geom_point(aes(fill = -log10(pvalue), size = gene_r), 
	       shape = 21, stroke = .5) +
    scale_y_reordered(position = "right") +
    scale_fill_gradient(low = "white", high = "black") +
    scale_size(range = c(.5, 3.5)) +
    facet_wrap(~module,
	       ncol = 1, scale = "free_y") +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  axis.text.x = element_blank(),
	  legend.text = element_text(size = 9, margin = margin(l = -0, unit = "lines")),
	  legend.title = element_text(size = 9),
	  legend.key.size = unit(.25, "cm"),
	  strip.text = element_blank(),
	  panel.grid = element_blank(),
	  plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")
	  ) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colorbar("logP:", barheight = 5, barwidth = .25),
	   size = guide_legend("Gene\nRatio:")) +
    coord_cartesian(clip = "off")

fig_c_title <- 
    ggdraw() + 
    draw_label(
	       "DN2 condition: gene programs are enriched in\ndifferent biological processes",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(-1, 1.25, -1, 2, unit = "lines"))

fig_c_grid <-
    plot_grid(
	      fig_c_title,
	      plot_grid(dn2_modules_plot, kim_plot, NULL, go_plot,
			nrow = 1, align = "h", rel_widths = c(.6, .45, .1, 1)),
	      ncol = 1, 
	      rel_heights = c(.1, 1),
	      labels = c(NULL, "c"), label_size = 10)

left_grid <-
    plot_grid(
	      fig_a_grid, NULL, fig_c_grid,
	      ncol = 1, rel_heights = c(0.5, .025, 1)
	      )


# Right-hand side ##############################################################

#### Test
stims <- c("Unstim", "IL4", "CD40L", "TLR9", "TLR7", "BCR", "BCR-TLR7", "DN2")

diff_expr <- 
    "./bcell_lowinput/results/edger/diff_expr_all_times.tsv" |>
    read_tsv()

diff_expr_summ <-
    diff_expr |>
    group_by(group1, group2) |>
    summarise(n = sum(!is.na(gene_id))) |>
    ungroup() |>
    separate(group1, c("stim1", "t1"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    separate(group2, c("stim2", "t2"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    mutate_at(vars(stim1, stim2), ~recode(., "BCR_TLR7" = "BCR-TLR7")) |>
    mutate_at(vars(stim1, stim2), ~factor(., levels = stims)) |>
    complete(group1, group2, fill = list(n = NA)) |>
    arrange(stim1, stim2, t1, t2) |>
    mutate_at(vars(group1, group2), ~factor(., levels = unique(c(group1, group2))))

x_lines_df <- 
    diff_expr_summ |>
    filter(!is.na(n)) |>
    distinct(group1, stim1, t1) |>
    rowid_to_column() |>
    group_by(stim1) |>
    summarise(avg = mean(rowid), 
	      rowid = max(rowid) + 0.5) |>
    ungroup()

y_lines_df <- 
    diff_expr_summ |>
    filter(!is.na(n)) |>
    distinct(group2, stim2, t2) |>
    rowid_to_column() |>
    group_by(stim2) |>
    summarise(avg = mean(rowid), 
	      rowid = max(rowid) + 0.5) |>
    ungroup()

testp <- 
    ggplot(data = diff_expr_summ, 
       aes(x = group1, y = group2)) +
    geom_tile(aes(fill = n), alpha = .9) +
    geom_vline(xintercept = head(x_lines_df$rowid, -1), 
	       color = "midnightblue", linewidth = .35) +
    geom_hline(yintercept = head(y_lines_df$rowid, -1), 
	       color = "midnightblue", linewidth = .35) +
    scale_x_discrete(labels = function(x) str_extract(x, "\\d+$")) +
    scale_y_discrete(labels = function(x) str_extract(x, "\\d+$")) +
    scale_fill_viridis_c(option = "cividis",
			 na.value = "white",
			 labels = scales::comma) +
    theme_minimal() +
    theme(
	  axis.title = element_blank(),
	  panel.grid.major = element_line(color = "black", linewidth = .35),
	  legend.position = "top",
	  legend.title.position = "top",
	  plot.margin = margin(b = 0.25, l = 0.75, unit = "in"),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(fill = "Number of DE genes:") +
    guides(fill = guide_colorbar(barwidth = 15, barheight = .5)) +
    annotate("text", x = x_lines_df$avg, y = -1, label = x_lines_df$stim1, 
	     size = 9, size.unit = "pt") +
    annotate("text", x = -1, y = y_lines_df$avg, label = y_lines_df$stim2, 
	     size = 9, size.unit = "pt", hjust = 1) +
    coord_cartesian(xlim = c(1, 28), ylim = c(1, 28), clip = "off")

ggsave("./paper_plots/testp.png", testp, width = 6, height = 6) 


########


# Fig B

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
    theme(text = element_text(size = 9),
	  axis.text = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  panel.grid.major.x = element_line(color = "grey90", linewidth = .2),
	  panel.grid.major.y = element_blank(),
	  strip.text.y = element_text(angle = 0, hjust = 0,
				      margin = margin(l = 0)
				      ),
	  panel.spacing = unit(0.1, "lines")
	  ) +
    labs(x = "hours", y = "Eigengene expression")

fig_b_title <- 
    ggdraw() + 
    draw_label(
	       "Gene programs in B cells",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 2, r = 0, unit = "lines"))

fig_b_grid <- 
    plot_grid(
	      fig_b_title, wgcna_plot,
	      ncol = 1, rel_heights = c(.1, 1),
	      labels = c(NULL, "b"), label_size = 10
    )


# Fig D SLE genes
sle_genes <- 
    c("PTPN22", "FCGR2A", "TNFSF4", "IL10", "NCF2", "SPRED2", "IFIH1", "STAT1", "STAT4",
      "IKZF1", "IKZF2", "IKZF3", "PXK", "IL12A", "BANK1", "BLK", "MIR146A", 
      "TNFAIP3", "IRF5", "IRF7", "IRF8", "TNIP1", "ATG5",
      "WDFY4", "ARID5B", "CD44", "ETS1", "SLC15A4", "CSK", "SOCS1", 
      "CLEC16A", "ITGAM", "TYK2", "UBE2L3", "TLR7", "IRAK1", "IKBKE")

gene_names <- 
    read_tsv("./bcell_lowinput/results/edger/results.tsv") |>
    mutate(stim = factor(stim, levels = rev(levels(legend_df$stim)))) |>
    distinct(gene_id, gene_name)

kme_df <- 
    "./bcell_lowinput/wgcna/data/DN2_kme.tsv" |>
    read_tsv() |>
    pivot_longer(-gene_id, names_to = "module", values_to = "kme") |>
    filter(module != "grey") |>
    left_join(gene_names, join_by(gene_id)) |>
    left_join(filter(module_sizes, stim == "DN2"), join_by(module)) |>
    select(gene_name, module = ix, kme) |>
    mutate(module = factor(module))

kme_sle_df <- 
    kme_df |>
    filter(gene_name %in% sle_genes) |>
    group_by(gene_name) |>
    filter(any(kme >= 0.8)) |>
    mutate(max_kme_mod = module[which.max(kme)]) |>
    ungroup() |>
    arrange(max_kme_mod, desc(kme)) |>
    mutate(gene_name = fct_inorder(gene_name))

sle_heatmap <- 
    ggplot(data = kme_sle_df, 
       aes(x = module, y = gene_name)) +
    geom_tile(aes(fill = kme)) +
    scale_fill_gradient2(low = "Navy Blue",
			 mid = "white",
			 high = "firebrick",
			 limits = c(-1, 1),
			 breaks = c(-1, 0, 1),
			 labels = c("-1", "0.5", "1")) +
    theme_minimal() +
    theme(text = element_text(size = 9)) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colorbar("kME:", barheight = 5, barwidth = .5))

fig_d_title <- 
    ggdraw() + 
    draw_label(
	       "Module membership of SLE\ngenes into 'DN2 modules'",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(0, 0, 0, 2, unit = "lines"))

fig_d_grid <- 
    plot_grid(fig_d_title, sle_heatmap,
	      ncol = 1, rel_heights = c(.1, 1),
	      labels = c("d", NULL), label_size = 10)


# Fig E LDSC
ldsc <- 
    "./bcell_lowinput/wgcna/gene_modules_LDSC.tsv" |>
    read_tsv() |>
    select(module = group, coefficient = Coefficient, std = Std, 
	   trait = Analysis, pfdr) |>
    mutate(trait = str_remove(trait, "_geneModules-DN2"),
	   trait = recode(trait, 
			  "AdultOnset" = "AO Asthma", 
			  "ChildOnset" = "CO Asthma",
			  "AllergyEczema" = "AllerEcz",
			  "Lupus" = "SLE",
			  "UKBAsthma" = "UO Asthma"),
	   trait = factor(trait, levels = c("Height", "UO Asthma", "CO Asthma", "AO Asthma", 
					    "AllerEcz", "PBC", "Arthritis", "SLE"))) |>
    left_join(filter(module_sizes, stim == "DN2"), join_by(module)) |>
    select(module = ix, trait, coefficient, std, pfdr) |>
    mutate(module = factor(module, levels = sort(unique(module))))


ldsc_plot <- 
    ggplot(data = ldsc,
	   aes(x = module, y = trait)) +
    geom_tile(aes(fill = coefficient)) +
    geom_point(data = filter(ldsc, pfdr <= 0.05),
	      aes(x = module, y = trait),
	      size = 1, shape = 8, color = "black") +
    scale_fill_gradient2(low = "Navy Blue", mid = "white", high = "firebrick") +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  axis.title = element_blank(),
	  panel.grid = element_blank()) +
    guides(fill = guide_colorbar("Coef:", barheight = 5, barwidth = .5))

fig_e_title <- 
    ggdraw() + 
    draw_label(
	       "DN2 module #3 is enriched in\nSLE heritability",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(0, 0, 0, 2, unit = "lines"))

fig_e_grid <- 
    plot_grid(fig_e_title, ldsc_plot,
	      ncol = 1, rel_heights = c(.1, 1),
	      labels = c("e", NULL), label_size = 10)

right_grid <- 
    plot_grid(
	      fig_b_grid, NULL, fig_d_grid, NULL, fig_e_grid, 
	      ncol = 1, rel_heights = c(.7, .05, 1, .05, 1)
	      )


ggsave("./paper_plots/fig2.png", 
       plot_grid(left_grid, NULL, right_grid, nrow = 1, rel_widths = c(1, .1, 0.625)),
       width = 6.5, height = 6.5, dpi = 600)

ggsave("./paper_plots/fig2.pdf", 
       plot_grid(left_grid, NULL, right_grid, nrow = 1, rel_widths = c(1, .1, 0.625)),
       width = 6.5, height = 6.5)


## Fig D Example Hub genes
#
#cpm_df <- 
#    "./bcell_lowinput/results/edger/cpm.tsv" |>
#    read_tsv() |>
#    group_by(sample_id, stim) |>
#    nest() |>
#    ungroup() |>
#    separate(sample_id, c("donor_id", "dummy", "time"), sep = "_") |>
#    mutate(hours = parse_number(time),
#	   hours = factor(hours, levels = sort(unique(hours))),
#           condition = paste(dummy, time, sep = " "),
#	   condition = factor(condition, levels = names(stim_colors)),
#	   stim = factor(stim, levels = rev(levels(legend_df$stim)))) |>
#    unnest(cols = data) |>
#    select(donor_id, condition, stim, hours, gene_id, gene_name, obs_cpm, obs_logcpm)
#
#
#genes_dn2_mod7 <- c("TK1", "DLGAP5", "CDCA5") 
#genes_tlr9_mod7 <- c("RPL36", "RPL32", "RPS3A")
#
#cpm_plot_df <-
#    cpm_df |>
#    filter(gene_name %in% c(genes_dn2_mod7, genes_tlr9_mod7),
#	   stim != "IL4")
#
#mod7_dn2 <- 
#    ggplot(data = cpm_plot_df |> 
#	   filter(gene_name %in% genes_dn2_mod7) |>
#	   mutate(gene_name = factor(gene_name, levels = genes_dn2_mod7)),
#	   aes(x = hours, y = obs_cpm)) +
#    geom_quasirandom(aes(color = condition),
#		     method = "smiley", width = .2, 
#		     size = .5, alpha = .7) +
#    geom_line(aes(group = stim), 
#	      stat = "smooth", method = "loess", span = 1, se = FALSE, 
#              linewidth = .7) +
#    scale_y_continuous(breaks = scales::pretty_breaks(2)) +
#    scale_color_manual(values = stim_colors) + 
#    facet_grid(gene_name ~ stim, scale = "free_y", drop = TRUE) +
#    theme_minimal() +
#    theme(text = element_text(size = 9),
#	  axis.text.x = element_blank(),
#	  panel.grid.minor = element_blank(),
#	  panel.grid.major.x = element_blank(),
#	  panel.grid.major.y = element_line(linetype = 3, 
#					    linewidth = .2,
#					    color = "grey66"),
#	  strip.text.x = element_blank(),
#	  strip.text.y = element_text(size = 8, 
#				      angle = 0, hjust = 1,
#				      margin = unit(c(0, 0, 0, 1), "pt")),
#	  plot.margin = unit(c(11, 5.5, 5.5, 5.5), "pt")
#	  ) +
#    labs(x = NULL, y = "CPM") +
#    guides(color = "none") +
#    coord_cartesian(clip = "off")
#
#mod7_tlr9 <- 
#    ggplot(data = cpm_plot_df |> 
#	   filter(gene_name %in% genes_tlr9_mod7) |>
#	   mutate(gene_name = factor(gene_name, levels = genes_tlr9_mod7)),
#	   aes(x = hours, y = obs_cpm)) +
#    geom_quasirandom(aes(color = condition),
#		     method = "smiley", width = .2, 
#		     size = .5, alpha = .7) +
#    geom_line(aes(group = stim), 
#	      stat = "smooth", method = "loess", span = 1, se = FALSE, 
#              linewidth = .7) +
#    scale_y_continuous(breaks = scales::pretty_breaks(2)) +
#    scale_color_manual(values = stim_colors) + 
#    facet_grid(gene_name ~ stim, scale = "free_y", drop = TRUE) +
#    theme_minimal() +
#    theme(text = element_text(size = 9),
#	  axis.text.x = element_blank(),
#	  panel.grid.minor = element_blank(),
#	  panel.grid.major.x = element_blank(),
#	  panel.grid.major.y = element_line(linetype = 3, 
#					    linewidth = .2,
#					    color = "grey66"),
#	  strip.text.x = element_blank(),
#	  strip.text.y = element_text(size = 8, 
#				      angle = 0, hjust = 1,
#				      margin = unit(c(0, 0, 0, 1), "pt")),
#	  plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")
#	  ) +
#    labs(x = NULL, y = "CPM") +
#    guides(color = "none") +
#    coord_cartesian(clip = "off")
#
#
#fig_d_title_1 <- 
#    ggdraw() + 
#    draw_label(
#	       "Genes in DN2 module #7",
#	       x = 0,
#	       size = 9,
#	       hjust = 0
#	       ) +
#    theme(text = element_text(size = 9),
#	  plot.margin = margin(1, 1.25, 1, 2, unit = "lines"))
#
#fig_d_title_2 <- 
#    ggdraw() + 
#    draw_label(
#	       "Genes in TLR9 module #7",
#	       x = 0,
#	       size = 9,
#	       hjust = 0
#	       ) +
#    theme(text = element_text(size = 9),
#	  plot.margin = margin(1, 1.25, 1, 2, unit = "lines"))
#
#fig_d_grid <- 
#    plot_grid(fig_d_title_1,
#	      mod7_dn2,
#	      fig_d_title_2,
#	      mod7_tlr9, 
#	      ncol = 1, rel_heights = c(.1, 1, .1, 1),
#	      labels = c("D)", NULL, NULL), label_size = 10)
#


























