library(tidyverse)
library(glue)
library(readxl)
library(cowplot)
#library(edgeR)
#library(fgsea)
library(mirai)

################################################################################
# Low-input data
################################################################################
stim_colors <- 
    "../paper_plots/figure_colors.txt" |>
    read_tsv(col_names = c("stim", "time", "color")) |>
    unite("condition", c(stim, time), sep = "_") |>
    mutate(condition = paste0(condition, "hrs")) |>
    deframe()

all_stims <-
    names(stim_colors) |>
    str_split("_") |>
    map_chr(1) |>
    unique()

diff_expr <- 
    "../bcell_lowinput/results/edger/diff_expr_all_times_all_genes.tsv" |>
    vroom::vroom() |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"),
	   stat = -log10(PValue) * sign(logFC)) |>
    arrange(contrast, stat)

vs_unstim0 <- 
    diff_expr |>
    filter(grepl("Unstim\\.0$", contrast)) |>
    separate(contrast, c("group1", "group2"), sep = "-") |>
    separate(group1, c("stim1", "t1"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    mutate(stim1 = recode(stim1, 
			  "IL4" = "IL-4c", 
			  "CD40L" = "CD40c", 
			  "TLR9" = "TLR9c", 
			  "TLR7" = "TLR7c",
			  "BCR" = "BCRc", 
			  "BCR_TLR7" = "BCR/TLR7c", 
			  "DN2" = "DN2c"),
	   stim1 = factor(stim1, levels = all_stims), 
	   condition = glue("{stim1}_{t1}hrs"),
	   condition = factor(condition, levels = names(stim_colors))) |>
    select(condition, stim1, t1, gene_id, gene_name, logFC, FDR) |>
    arrange(condition, FDR)

################################################################################
# Perez
################################################################################
perez <- 
    read_tsv("./data/science.abf1970_table_s7", 
	     skip = 1, 
	     col_names = c("cell", "gene", "p_val", "beta", "se")) |>
    filter(cell == "b")

# Yazar
yazar <- 
    read_excel("./data/science.abf3041_tables_s6_to_s19.xlsx", 
	       sheet = "Table.S10", 
	       skip = 2) |>
    janitor::clean_names() |>
    filter(cell_type %in% c("B IN", "B Mem")) |>
    filter(e_snp_rank == "eSNP1") |>
    group_by(gene_id, gene_ensembl_id) |>
    slice_min(pvalue, n = 1, with_ties = FALSE) |>
    ungroup()

perez_genes <- 
    anti_join(perez, yazar, join_by(gene == gene_id)) |>
    filter(p_val < 5e-5) |>
    mutate(gene = recode(gene, 
			 "HIST1H2AE" = "H2AC8",
			 "ELMSAN1" = "MIDEAS",
			 "HIST1H4C" = "H4C3",
			 "SELT" = "SELENOT",
			 "ZNRD1" = "POLR1H",
			 "TMEM66" = "SARAF",
			 "PCNX" = "PCNX1",
			 "METTL10" = "EEF1AKMT2",
			 "TRAPPC2P1" = "TRAPPC2",
			 "FAM63B" = "MINDY2",
			 "C5orf45" = "MRNIP"))

# Proportion of genes DE in our low-input data
#perez <- 
#    "./data/Bcell_eGenes_inPerezButNotYazar_so_SLEenriched.txt" |>
#    read_tsv(col_names = "gene_name") 
    
vs_unstim0_perez_gene <-
    perez_genes |>
    select(gene_name = gene) |>
    expand_grid(distinct(vs_unstim0, condition, stim1, t1)) |>
    left_join(vs_unstim0, join_by(gene_name, condition, stim1, t1)) |>
    mutate(category = case_when(logFC > 0 & FDR <= 0.05 ~ "up",
				logFC < 0 & FDR <= 0.05 ~ "down",
				FDR > 0.05 ~ "n.s.",
				is.na(FDR) ~ "N.A.",
				TRUE ~ NA_character_))

vs_unstim0_perez_data <-
    vs_unstim0_perez_gene |>
    count(condition, stim1, t1, category) |>
    mutate(category = fct_inorder(category))

vs_unstim0_perez_plot <- 
    ggplot(vs_unstim0_perez_data) +
    geom_col(aes(x = factor(t1), y = n, fill = category), 
             color = "black", position = "fill", linewidth = .25) +
    geom_hline(yintercept = .5, linetype = 2, linewidth = 0.25, color = "black") +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    scale_fill_manual(values = c("grey", "001260", "white", "#AA4613")) +
    facet_grid(.~stim1, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          strip.text.x = element_text(size = 7),
          strip.text.y = element_text(size = 7),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9)
    ) +
    labs(x = "Hours", y = "Proportion of gene signature", fill = "DE direction:",
	 title = "Perez et al. B cell eGenes")

ggsave("vs_unstim_perez.png",
       vs_unstim0_perez_plot, 
       width = 8, height = 2)

vs_unstim0_perez_data |>
    group_by(condition) |>
    mutate(p = n/sum(n) * 100) |>
    ungroup() |>
    filter(category == "up") |>
    group_by(stim1) |>
    summarise(m = max(p)) |>
    ungroup()

# upset plot
# perez_genes_up_list <- 
#     vs_unstim0_perez_gene |>
#     filter(category == "up") |>
#     select(condition, gene_name) |>
#     {function(x) split(x, x$condition)}() |>
#     map(~pull(., gene_name))
# 
# library(UpSetR)
# 
# ?upset
# 
# png("./upset.png", width = 7, height = 7, units = "in", res = 300)
# upset(fromList(perez_genes_up_list), 
#       nsets = length(perez_genes_up_list), 
#       keep.order = TRUE)
# dev.off()
# 
# 
# vs_unstim0_perez_gene |>
#     filter(category == "up") |>
#     select(condition, gene_name) |>
#     count(gene_name) |>
#     arrange(desc(n)) |>
#     print(n = Inf)




# Enrichment
perez_genes_ids <- 
    perez_genes |> 
    inner_join(gencode_v38, join_by(gene == gene_name))	

perez_list <- list(perez = perez_genes_ids$gene_id)

diff_expr <- 
    "../bcell_lowinput/results/edger/diff_expr_all_times_all_genes.tsv" |>
    vroom::vroom() |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"),
	   stat = -log10(PValue) * sign(logFC)) |>
    arrange(contrast, stat)

bcell_gene_lists <-
    diff_expr |>
    select(contrast, SYMBOL = gene_id, stat) |>
    {function(x) split(x, x$contrast)}() |>
    map(~select(., -contrast) |> deframe())

daemons(16)

gsea_res_perez <- 
    map(bcell_gene_lists,
	in_parallel(function(x) fgsea::fgsea(pathways = gene_sets, stat = x, nPermSimple = 100000), 
		    gene_sets = perez_list))

daemons(0)

gsea_res_perez_df <- 
    gsea_res_perez |>
    bind_rows(.id = "contrast") |>
    as_tibble()

write_tsv(gsea_res_perez_df, "./data/gsea_results_perez.tsv")
write_rds(gsea_res_perez, "./data/gsea_results_perez.rds")

# plot
stim_colors <- 
    "../paper_plots/figure_colors.txt" |>
    read_tsv(col_names = c("stim", "time", "color")) |>
    unite("condition", c(stim, time), sep = "_") |>
    mutate(condition = paste0(condition, "hrs")) |>
    deframe()

all_stims <-
    names(stim_colors) |>
    str_split("_") |>
    map_chr(1) |>
    unique()

diff_expr_summ <-
    gsea_res_perez_df |>
    mutate(padj_global = p.adjust(pval, method = "fdr")) |>
    separate(contrast, c("group1", "group2"), sep = "-") |>
    separate(group1, c("stim1", "t1"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    separate(group2, c("stim2", "t2"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    mutate_at(vars(stim1, stim2), ~recode(., 
					 "IL4" = "IL-4c", 
					 "CD40L" = "CD40c", 
					 "TLR9" = "TLR9c", 
					 "TLR7" = "TLR7c",
					 "BCR" = "BCRc", 
					 "BCR_TLR7" = "BCR/TLR7c", 
					 "DN2" = "DN2c")) |>
    mutate_at(vars(stim1, stim2), ~factor(., levels = all_stims)) |>
    select(group1, stim1, t1, group2, stim2, t2, pathway, pval, padj = padj_global, NES) |>
    arrange(stim1, stim2, t1, t2) |>
    mutate_at(vars(group1, group2), ~factor(., levels = unique(c(group1, group2))))

segments_y <- 
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

segments_x <- 
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

b_axis_y_df <- 
    diff_expr_summ |>
    filter(!is.na(stim1)) |>
    distinct(group1, stim1, t1) |>
    mutate(condition1 = glue("{stim1}_{t1}hrs")) |>
    select(group1, condition1)

b_axis_x_df <- 
    diff_expr_summ |>
    filter(!is.na(stim2)) |>
    distinct(group2, stim2, t2) |>
    mutate(condition2 = glue("{stim2}_{t2}hrs")) |>
    select(group2, condition2)

b_axis_y_plot <-
    ggplot(data = b_axis_y_df, 
           aes(x = factor(1), y = group1)) +
    geom_tile(aes(fill = condition1)) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(), 
          plot.margin = margin(0, 0, 0, 0)
    ) +
    guides(fill = "none") +
    coord_cartesian(ylim = c(1, 28), xlim = c(1, 1))

b_axis_x_plot <-
    ggplot(data = b_axis_x_df, 
           aes(y = factor(1), x = group2)) +
    geom_tile(aes(fill = condition2)) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(), 
          plot.margin = margin(0, 0, 0, 0)
    ) +
    guides(fill = "none") +
    coord_cartesian(ylim = c(1, 1), xlim = c(1, 28))

diff_plot <- 
    ggplot(data = diff_expr_summ) +
    geom_tile(aes(x = group2, y = group1, fill = NES)) +
    geom_vline(xintercept = segments_x$rowid, color = "white", linewidth = .5) +
    geom_hline(yintercept = segments_y$rowid, color = "white", linewidth = .5) +
    geom_text(data = filter(diff_expr_summ, padj < 0.05),
	      aes(x = group2, group1, label = "*"), 
	      color = "white", alpha = .75, fontface = "bold") +
    scico::scale_fill_scico(palette = "vik", 
			    na.value = "white",
			    limits = c(-3, 3),
			    labels = scales::comma) +
    theme_minimal() +
    theme(text = element_text(size = 9),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  panel.grid.major = element_blank(),
	  panel.background = element_rect(color = NA, fill = "transparent"),
	  plot.margin = margin(0, 0, 0, 0),
	  legend.position.inside = c(.8, .3),
	  ) +
    labs(fill = "NES:") +
    guides(fill = guide_colorbar(barwidth = .5, barheight = 5, position = "inside")) +
    coord_cartesian(xlim = c(1, 28), ylim = c(1, 28))

plot_tmp <- 
    plot_grid(NULL, b_axis_x_plot, b_axis_y_plot, diff_plot,
	      ncol = 2, nrow = 2, align = "v", 
	      rel_heights = c(.075, 1), rel_widths = c(0.06, 1))

enrich_plot <- 
    plot_grid(
	      plot_title,
	      plot_tmp + theme(plot.margin = margin(t = .5, unit = "lines")),
	      ncol = 1, 
	      rel_heights = c(.1, 1)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./enrichments_perez.png", enrich_plot, width = 4, height = 4)




################################################################################
# Gencode annotations
################################################################################
#gencode_v38 <- 
#    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
#              "gencode.v38.primary_assembly.annotation.gtf.gz") |>
#    read_tsv(comment = "#", col_names = FALSE) |>
#    filter(X3 == "gene") |>
#    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
#           gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
#    select(gene_id, gene_name) |>
#    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))


# In Nakano et al, they say they followed analysis from Ota et al 2021, where they used Gencode v27
# Still, some gene names don't match those annotated in Gencode v27
gencode_v27 <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
              "gencode.v27.primary_assembly.annotation.gtf.gz") |>
    rtracklayer::import(feature.type = "gene") |>
    as.data.frame() |>
    as_tibble() |>
    select(gene_id, gene_name) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

dis_state <- 
    read_excel("./data/mmc2.xlsx", sheet = 2, skip = 2) |>
    janitor::clean_names()

dis_activ <- 
    read_excel("./data/mmc2.xlsx", sheet = 3, skip = 2) |>
    janitor::clean_names()

# Proportion of gene signature captured by up-regulated genes in our data
set_state_b <- 
    dis_state |>
    filter(cell_type %in% c("Naive B", "USM B", "SM B", "DN B", "Plasmablast"),
	   !grepl("^\\d+$", gene),
	   log_fc > 0, fdr <= 0.05) |>
    mutate(cell_type = factor(cell_type, levels = c("Naive B", "USM B", "SM B", "DN B", "Plasmablast"))) |>
    left_join(gencode_v27, join_by(gene == gene_name), relationship = "many-to-many") |>
    group_by(gene) |>
    mutate(n = n_distinct(gene_id)) |>
    ungroup() |>
    filter(!is.na(gene_id)) |>
    select(cell_type, gene_id)

set_activ_b <- 
    dis_activ |>
    filter(cell_type %in% c("Naive B", "USM B", "SM B", "DN B", "Plasmablast"),
	   !grepl("^\\d+$", gene),
	   log_fc > 0, fdr <= 0.05) |>
    mutate(cell_type = factor(cell_type, levels = c("Naive B", "USM B", "SM B", "DN B", "Plasmablast"))) |>
    left_join(gencode_v27, join_by(gene == gene_name), relationship = "many-to-many") |>
    group_by(gene) |>
    mutate(n = n_distinct(gene_id)) |>
    ungroup() |>
    filter(!is.na(gene_id)) |>
    select(cell_type, gene_id)

vs_unstim0_state_data <- 
    set_state_b |>
    expand_grid(distinct(vs_unstim0, condition, stim1, t1)) |>
    left_join(vs_unstim0, join_by(gene_id, condition, stim1, t1)) |>
    mutate(category = case_when(logFC > 0 & FDR <= 0.05 ~ "up",
				logFC < 0 & FDR <= 0.05 ~ "down",
				FDR > 0.05 ~ "n.s.",
				is.na(FDR) ~ "N.A.",
				TRUE ~ NA_character_)) |>
    count(cell_type, condition, stim1, t1, category) |>
    mutate(category = fct_inorder(category))

vs_unstim0_activ_data <- 
    set_activ_b |>
    expand_grid(distinct(vs_unstim0, condition, stim1, t1)) |>
    left_join(vs_unstim0, join_by(gene_id, condition, stim1, t1)) |>
    mutate(category = case_when(logFC > 0 & FDR <= 0.05 ~ "up",
				logFC < 0 & FDR <= 0.05 ~ "down",
				FDR > 0.05 ~ "n.s.",
				is.na(FDR) ~ "N.A.",
				TRUE ~ NA_character_)) |>
    count(cell_type, condition, stim1, t1, category) |>
    mutate(category = fct_inorder(category))

nakano_vs0_data <-
    bind_rows("state" = vs_unstim0_state_data, "activity" = vs_unstim0_activ_data, .id = "set") |>
    mutate(set = fct_inorder(set)) |>
    arrange(set, cell_type, stim1, t1, category) |>
    unite("combined_group", c(set, category), sep = "_", remove = FALSE) |>
    mutate(combined_group = fct_inorder(combined_group),
	   x_dodged = ifelse(set == "state", as.numeric(factor(t1)) - 0.17, as.numeric(factor(t1)) + 0.17))

vs0_colors <- 
    c(
      "state_up" = "#800020",
      "activity_up" = "#CD5C5C",
      "state_down" = "royalblue",
      "activity_down" = "#86BBD8",
      "state_n.s." = "gray50",
      "activity_n.s." = "gray70",
      "state_N.A." = "white",
      "activity_N.A." = "white")

legend_colors <- c("up" = "#800020", "down" = "royalblue", "n.s." = "grey50", "N.A.")
legend_shades <- c("state" = "#beb3a4", "activity" = "#e9d4b7")


nakano_plot <- 
    ggplot(nakano_vs0_data, aes(x = x_dodged, y = n)) +
    geom_col(aes(fill = combined_group), 
             color = "black", position = "fill", linewidth = .2, width = .3) +
    scale_fill_manual(values = vs0_colors) +
    scale_x_continuous(
		       breaks = 1:4, 
		       labels = as.character(unique(nakano_vs0_data$t1))) +
    scale_y_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
    geom_point(aes(y = 0.5, color = category), alpha = 0) +
    scale_color_manual(name = "Regulation", values = legend_colors) +
    geom_point(aes(y = 0.5, shape = set), alpha = 0) +
    scale_shape_manual(name = "Gene Set", values = c(15, 15)) +
    facet_grid(cell_type~stim1, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(axis.text = element_text(size = 7),
          axis.title = element_text(size = 8),
          strip.text.x = element_text(size = 6),
          strip.text.y = element_text(size = 6),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing = unit(0.2, "lines"),
          legend.box.margin = margin(l = -5),
          legend.margin = margin(0, 0, 0, 0),
          strip.background = element_rect(linewidth = 0.2, color = "black"),
          panel.border = element_rect(linewidth = 0.2, fill = NA),
    ) +
    labs(x = "Hours", y = "Proportion of gene signature") +
    guides(
        fill = "none", 
        color = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 5), order = 1),
        shape = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 5, color = legend_shades), order = 2)
    )

ggsave("nakano_signatures.png", nakano_plot, width = 6.5, height = 3.5)


nakano_vs0_data |>
    group_by(set, cell_type, stim = stim1, timep = t1) |>
    mutate(p = n/sum(n) * 100) |>
    ungroup() |>
    filter(set == "activity", category == "up") |>
    select(cell_type, stim, timep, p) |>
    group_by(stim) |>
    summarise(min = min(p), max = max(p)) |>
    ungroup() |>
    filter(stim != "Unstim") |>
    summarise(min(min), max(max))

nakano_vs0_data |>
    group_by(set, cell_type, stim = stim1, timep = t1) |>
    mutate(p = n/sum(n) * 100) |>
    ungroup() |>
    filter(set == "activity", category == "up") |>
    select(cell_type, stim, timep, p) |>
    filter(cell_type == "Plasmablast") |>
    group_by(stim) |>
    summarise(min = min(p), max = max(p)) |>
    ungroup()




# separate plots for each gene set:
vs_unstim0_state_plot <- 
    ggplot(vs_unstim0_state_data) +
    geom_col(aes(x = factor(t1), y = n, fill = category), 
             color = "black", position = "fill", linewidth = .25) +
    geom_hline(yintercept = .5, linetype = 2, linewidth = 0.25, color = "black") +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    scale_fill_manual(values = c("grey", "001260", "white", "#AA4613")) +
    facet_grid(cell_type~stim1, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          strip.text.x = element_text(size = 7),
          strip.text.y = element_text(size = 7),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9)
    ) +
    labs(x = "Hours", y = "Proportion of gene signature", fill = "DE direction:",
	 title = "Disease-state signatures")

vs_unstim0_activ_plot <- 
    ggplot(vs_unstim0_activ_data) +
    geom_col(aes(x = factor(t1), y = n, fill = category), 
             color = "black", position = "fill", linewidth = .25) +
    geom_hline(yintercept = .5, linetype = 2, linewidth = 0.25, color = "black") +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    scale_fill_manual(values = c("grey", "001260", "white", "#AA4613")) +
    facet_grid(cell_type~stim1, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          strip.text.x = element_text(size = 7),
          strip.text.y = element_text(size = 7),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9)
    ) +
    labs(x = "Hours", y = "Proportion of gene signature", fill = "DE direction:",
	 title = "Disease-activity signatures")

ggsave("vs_unstim.png",
       plot_grid(vs_unstim0_state_plot, vs_unstim0_activ_plot, ncol = 1), 
       width = 8, height = 8)








# Enrichment
gene_sets <-
    dis_state |> 
    filter(cell_type %in% c("Naive B", "USM B", "SM B", "DN B", "Plasmablast"),
	   log_fc > 0, fdr < 0.05) |>
    left_join(gencode_v27, join_by(gene == gene_name)) |>
    group_by(gene) |>
    mutate(n = n_distinct(gene_id)) |>
    ungroup() |>
    filter(!is.na(gene_id)) |>
    {function(x) split(x, x$cell_type)}() |>
    map(~arrange(., log_fc) |> pull(gene_id))

diff_expr <- 
    "../bcell_lowinput/results/edger/diff_expr_all_times_all_genes.tsv" |>
    vroom::vroom() |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"),
	   stat = -log10(PValue) * sign(logFC)) |>
    arrange(contrast, stat)

bcell_gene_lists <-
    diff_expr |>
    select(contrast, SYMBOL = gene_id, stat) |>
    {function(x) split(x, x$contrast)}() |>
    map(~select(., -contrast) |> deframe())

daemons(16)

gsea_res <- 
    map(bcell_gene_lists,
	in_parallel(function(x) fgsea::fgsea(pathways = gene_sets, stat = x, nPermSimple = 100000), 
		    gene_sets = gene_sets))

daemons(0)

gsea_res_df <- 
    gsea_res |>
    bind_rows(.id = "contrast") |>
    as_tibble()

write_tsv(gsea_res_df, "./data/gsea_results.tsv")
write_rds(gsea_res, "./data/gsea_results.rds")



# Plot
gsea_res_df <- read_tsv("./data/gsea_results.tsv")
gsea_res <- read_rds("./data/gsea_results.rds")

stim_colors <- 
    "../paper_plots/figure_colors.txt" |>
    read_tsv(col_names = c("stim", "time", "color")) |>
    unite("condition", c(stim, time), sep = "_") |>
    mutate(condition = paste0(condition, "hrs")) |>
    deframe()

all_stims <-
    names(stim_colors) |>
    str_split("_") |>
    map_chr(1) |>
    unique()

diff_expr_summ <-
    gsea_res_df |>
    mutate(padj_global = p.adjust(pval, method = "fdr")) |>
    separate(contrast, c("group1", "group2"), sep = "-") |>
    separate(group1, c("stim1", "t1"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    separate(group2, c("stim2", "t2"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    mutate_at(vars(stim1, stim2), ~recode(., 
					 "IL4" = "IL-4c", 
					 "CD40L" = "CD40c", 
					 "TLR9" = "TLR9c", 
					 "TLR7" = "TLR7c",
					 "BCR" = "BCRc", 
					 "BCR_TLR7" = "BCR/TLR7c", 
					 "DN2" = "DN2c")) |>
    mutate_at(vars(stim1, stim2), ~factor(., levels = all_stims)) |>
    select(group1, stim1, t1, group2, stim2, t2, pathway, pval, padj = padj_global, NES) |>
    arrange(stim1, stim2, t1, t2) |>
    mutate_at(vars(group1, group2), ~factor(., levels = unique(c(group1, group2))))

segments_y <- 
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

segments_x <- 
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

b_axis_y_df <- 
    diff_expr_summ |>
    filter(!is.na(stim1)) |>
    distinct(group1, stim1, t1) |>
    mutate(condition1 = glue("{stim1}_{t1}hrs")) |>
    select(group1, condition1)

b_axis_x_df <- 
    diff_expr_summ |>
    filter(!is.na(stim2)) |>
    distinct(group2, stim2, t2) |>
    mutate(condition2 = glue("{stim2}_{t2}hrs")) |>
    select(group2, condition2)

b_axis_y_plot <-
    ggplot(data = b_axis_y_df, 
           aes(x = factor(1), y = group1)) +
    geom_tile(aes(fill = condition1)) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(), 
          plot.margin = margin(0, 0, 0, 0)
    ) +
    guides(fill = "none") +
    coord_cartesian(ylim = c(1, 28), xlim = c(1, 1))

b_axis_x_plot <-
    ggplot(data = b_axis_x_df, 
           aes(y = factor(1), x = group2)) +
    geom_tile(aes(fill = condition2)) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(), 
          plot.margin = margin(0, 0, 0, 0)
    ) +
    guides(fill = "none") +
    coord_cartesian(ylim = c(1, 1), xlim = c(1, 28))


diff_plots <- 
    diff_expr_summ |>
    {function(x) split(x, x$pathway)}() |>
    map(function(plot_data) {
        
        diff_plot <- 
            ggplot(data = plot_data) +
            geom_tile(aes(x = group2, y = group1, fill = NES), alpha = .85) +
            geom_vline(xintercept = segments_x$rowid, color = "white", linewidth = .5) +
            geom_hline(yintercept = segments_y$rowid, color = "white", linewidth = .5) +
            geom_text(data = filter(plot_data, padj < 0.05),
                      aes(x = group2, group1, label = "*"), 
                      color = "white", alpha = .75, fontface = "bold") +
            scico::scale_fill_scico(palette = "vik", 
                                    na.value = "white",
                                    limits = c(-3, 3),
                                    labels = scales::comma) +
            theme_minimal() +
            theme(text = element_text(size = 9),
        	  axis.text = element_blank(),
        	  axis.title = element_blank(),
        	  panel.grid.major = element_blank(),
        	  panel.background = element_rect(color = NA, fill = "transparent"),
        	  plot.margin = margin(0, 0, 0, 0),
        	  legend.position.inside = c(.8, .3),
        	  ) +
            labs(fill = "NES:") +
            guides(fill = guide_colorbar(barwidth = .5, barheight = 5, position = "inside")) +
            coord_cartesian(xlim = c(1, 28), ylim = c(1, 28))


        plot_title <- 
            ggdraw() + 
            draw_label(
                unique(plot_data$pathway),
                x = 0,
                size = 9,
                hjust = 0
            ) +
            theme(text = element_text(size = 9),
                  plot.margin = margin(l = 1.25, unit = "lines"))

        plot_tmp <- 
            plot_grid(NULL, b_axis_x_plot, b_axis_y_plot, diff_plot,
                      ncol = 2, nrow = 2, align = "v", 
                      rel_heights = c(.075, 1), rel_widths = c(0.06, 1))

        plot_grid(
            plot_title,
            plot_tmp + theme(plot.margin = margin(t = .5, unit = "lines")),
            ncol = 1, 
            rel_heights = c(.1, 1)) +
        theme(plot.background = element_rect(fill = "white", color = "white"))
    })

plot_out <- plot_grid(plotlist = diff_plots, nrow = 2, ncol = 2)

ggsave("./enrichments.png", plot_out, width = 8, height = 8)

# Try limma camera
library(edgeR)

dge <- 
    read_rds("../bcell_lowinput/results/edger/diff_expr_all_times_dge.rds")

design <- 
    read_rds("../bcell_lowinput/results/edger/diff_expr_all_times_design.rds")

v <- voom(dge, design)
expr_mat <- v$E
rownames(expr_mat) <- str_remove(rownames(expr_mat), "\\.\\d+$")

idx_list <- 
    ids2indices(gene_sets, rownames(expr_mat), remove.empty = TRUE)

contrast_names <- combn(colnames(design), 2, simplify = FALSE)

contrast_matrix <- 
    makeContrasts(contrasts = sapply(contrast_names, function(x) paste(rev(x), collapse = "-")),
                  levels = design)

# camera_results <- 
#     lapply(seq_len(ncol(contrast_matrix)), 
#            function(i) {
#                cam_res <- camera(expr_mat, 
#                                  index = idx_list, 
#                                  design = design, 
#                                  contrast = contrast_matrix[ , i])
#                cam_res$contrast <- colnames(contrast_matrix)[i]
#                cam_res
#            })
# 
# camera_res_df <- 
#     camera_results |>
#     map_dfr(~rownames_to_column(., "cell_type") |> as_tibble())
# 
# write_tsv(camera_res_df, "./data/camera_results.tsv")

camera_res_df <- read_tsv("./data/camera_results.tsv")

camera_summ <-
    camera_res_df |>
    separate(contrast, c("group1", "group2"), sep = "-") |>
    separate(group1, c("stim1", "t1"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    separate(group2, c("stim2", "t2"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    mutate_at(vars(stim1, stim2), ~recode(., 
                                          "IL4" = "IL-4c", 
                                          "CD40L" = "CD40c", 
                                          "TLR9" = "TLR9c", 
                                          "TLR7" = "TLR7c",
                                          "BCR" = "BCRc", 
                                          "BCR_TLR7" = "BCR/TLR7c", 
                                          "DN2" = "DN2c")) |>
    mutate_at(vars(stim1, stim2), ~factor(., levels = all_stims)) |>
    select(group1, stim1, t1, group2, stim2, t2, cell_type, PValue, FDR, Direction) |>
    arrange(stim1, stim2, t1, t2) |>
    mutate_at(vars(group1, group2), ~factor(., levels = unique(c(group1, group2))))


diff_plots_camera <- 
    camera_summ |>
    {function(x) split(x, x$cell_type)}() |>
    map(function(plot_data) {
        
        diff_plot <- 
            ggplot(data = plot_data) +
            geom_tile(aes(x = group2, y = group1, fill = Direction)) +
            geom_vline(xintercept = segments_x$rowid, color = "white", linewidth = .5) +
            geom_hline(yintercept = segments_y$rowid, color = "white", linewidth = .5) +
            geom_text(data = filter(plot_data, FDR < 0.05),
                      aes(x = group2, group1, label = "*"), 
                      color = "white", alpha = .75, fontface = "bold") +
            scale_fill_manual(values = c("Down" = "001260", "Up" = "#AA4613"),
                              na.value = "white") +
            theme_minimal() +
            theme(text = element_text(size = 9),
                  axis.text = element_blank(),
                  axis.title = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.background = element_rect(color = NA, fill = "transparent"),
                  plot.margin = margin(0, 0, 0, 0),
                  legend.position.inside = c(.8, .3),
            ) +
            labs(fill = "Direction:") +
            coord_cartesian(xlim = c(1, 28), ylim = c(1, 28)) +
            guides(fill = guide_legend(position = "inside"))
        
        
        plot_title <- 
            ggdraw() + 
            draw_label(
                unique(plot_data$cell_type),
                x = 0,
                size = 9,
                hjust = 0
            ) +
            theme(text = element_text(size = 9),
                  plot.margin = margin(l = 1.25, unit = "lines"))
        
        plot_tmp <- 
            plot_grid(NULL, b_axis_x_plot, b_axis_y_plot, diff_plot,
                      ncol = 2, nrow = 2, align = "v", 
                      rel_heights = c(.075, 1), rel_widths = c(0.06, 1))
        
        plot_grid(
            plot_title,
            plot_tmp + theme(plot.margin = margin(t = .5, unit = "lines")),
            ncol = 1, 
            rel_heights = c(.1, 1)) +
            theme(plot.background = element_rect(fill = "white", color = "white"))
    })

plot_out_camera <- plot_grid(plotlist = diff_plots_camera, nrow = 2, ncol = 2)



ggsave("./enrichments_camera.png", plot_out_camera, width = 8, height = 8)


################################################################################
# TFs with motif enrichment by HOMER
################################################################################
#stims <- read_lines("../atacseq/homer/stims.txt")
#
#homer_df <-
#    glue("../atacseq/homer/results/{stims}/knownResults.txt") |>
#    setNames(paste0(stims, "c")) |>
#    map_dfr(~read_tsv(.) |>
#            janitor::clean_names() |>
#            {function(x) setNames(x, str_remove(names(x), "_of_\\d+$"))}(),
#            .id = "stim") |>
#    mutate_at(vars(starts_with("percent_of_")), parse_number) |>
#    mutate(stim = factor(stim, levels = paste0(stims, "c")),
#           fc = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) |>
#    extract(motif_name,
#	    c("tf_name", "tf_family", "dataset", "db"), 
#	    "(.+?)(?:\\((.+)\\))?/([^/]+)/(.+)", 
#	    remove = FALSE) |>
#    mutate(log10p = log_p_value/log(10)) |>
#    select(stim, motif_name, tf_name, tf_family, dataset, consensus,
#           pct_target =  percent_of_target_sequences_with_motif,
#           pct_bg = percent_of_background_sequences_with_motif,
#           fc, log10p, q_value = q_value_benjamini)
#
#homer_motifs <- 
#    homer_df |>
#    filter(q_value <= 0.01) |>
#    distinct(motif_name, tf_name, tf_family)
#
#write_tsv(homer_motifs, "./data/tf_motifs_fdr01.tsv")

# Plot
#tf_df <- 
#    "https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.txt" |>
#    read_tsv() |>
#    janitor::clean_names() 
#
#tf_enrich <- 
#    c("BATF", "ATF3", "FOS", "JUNB",  "FOSL1", "FOSL2", "BACH2", "NFE2L2", "BACH1", "NFE2",
#      "EBF1", "PBX2", "IRF1", "IRF2", "IRF4", "IRF8", "RELA", "RELB", "NFATC1", "NFATC2",
#      "RUNX1", "RUNX2", "STAT6", "STAT5A", "STAT5B", "STAT1", "TBX21", "PRDM1",
#      "EGR1", "EGR2", "BCL6")

tf_enrich <- 
    c("BATF", "ATF3", "FOSL1", "JUN", "JUNB", "FOSL2", "BACH2", "MAFK", "BACH1", "NFE2L2",
      "STAT6", "STAT5A", "STAT5B", "STAT1", "STAT4", "STAT3", 
      "REL", "RELA", "NFKB1", "NFKB2", "NFATC1", "NFATC2", 
      "PRDM1", "EGR1", "EGR2", "BCL6", "ZBTB33",
      "RUNX1", 
      "EBF1",
      "TBX21",
      "NR4A1",
      "SPI1",
      "POU2F2", "POU5F1",
      "ZNF143", "RBPJ",
      "MEF2A", "MEF2C")

cell_types_b <- c("Naive B", "USM B", "SM B", "DN B", "Plasmablast")

tf_nakano_data <- 
    bind_rows("Disease state" = dis_state, "Disease activity" = dis_activ, .id = "contrast") |> 
    filter(cell_type %in% cell_types_b, gene %in% tf_enrich) |>
    mutate(cell_type = factor(cell_type, levels = cell_types_b),
	   gene = factor(gene, levels = rev(tf_enrich)))

tf_state_plot <- 
    ggplot(tf_nakano_data) +
    geom_tile(aes(x = cell_type, y = gene, fill = log_fc)) +
    geom_text(data = filter(tf_nakano_data, fdr <= 0.05),
	      aes(x = cell_type, y = gene, label = "*"), 
	      vjust = .75) +
    scale_fill_gradient2(
			 low = "#2166AC",    # blue
			 mid = "white", 
			 high = "#B2182B",   # red
			 midpoint = 0,
			 limits = c(-2.5, 2.5),
			 oob = scales::squish
    ) +
    facet_wrap(~fct_inorder(contrast), nrow = 1) +
    theme_minimal() +
    theme(axis.text = element_text(size = 7),
	  axis.title = element_text(size = 8),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = "Cell type in Nakano et al.", y = NULL) +
    guides(fill = guide_colorbar(barwidth = .5))

ggsave("./enrich_tf_state.png", tf_state_plot, width = 6.5, height = 4)




# Manually curated to add gene and species
homer_annot_df <- read_excel("./data/tf_motifs_homer.xlsx")
    
homer_annot_df |>
    filter(species == "Homo sapiens") |>
    select(motif_name, Gene) |>
    separate_rows(Gene, sep = ",") |>
    left_join(gencode, join_by(Gene == gene_name)) |>
    filter(is.na(gene_id))

mouse_genes <- 
    homer_annot_df |>
    filter(species == "Mus musculus") |>
    separate_rows(Gene, sep = ",") |>
    drop_na() |>
    distinct(Gene)

library(babelgene)

# Get more detailed info including Ensembl IDs
library(babelgene)

orthologs_res <- 
    orthologs(
	      genes = mouse_genes$Gene,
	      species = "Mus musculus", 
	      human = FALSE,
	      top = FALSE
    ) |>
    as_tibble() |>
    distinct(human_symbol, human_ensembl, symbol, ensembl)

mouse_genes_orthos <- 
    left_join(mouse_genes, orthologs_res, join_by(Gene == symbol)) |>
    mutate(human_symbol = case_when(is.na(human_symbol) ~ toupper(Gene),
				    .default = human_symbol),
	   human_symbol = ifelse(human_symbol == "STAT5", "STAT5A", human_symbol))

human_genes <- 
    homer_annot_df |>
    filter(species == "Homo sapiens") |>
    select(Gene) |>
    separate_rows(Gene, sep = ",")

tfs_curated <- 
    mouse_genes_orthos |>
    select(Gene = human_symbol) |>
    bind_rows(human_genes)

# Check if these TFs are DE between patients and controls
dis_activ_b_all <- 
    read_excel("./data/mmc2.xlsx", sheet = 3, skip = 2) |>
    janitor::clean_names() |>
    filter(cell_type %in% c("Naive B", "USM B", "SM B", "DN B", "Plasmablast"),
	   !grepl("^\\d+$", gene)) |>
    mutate(cell_type = factor(cell_type, levels = c("Naive B", "USM B", "SM B", "DN B", "Plasmablast")))

dis_state_b_all <- 
    read_excel("./data/mmc2.xlsx", sheet = 2, skip = 2) |>
    janitor::clean_names() |>
    filter(cell_type %in% c("Naive B", "USM B", "SM B", "DN B", "Plasmablast"),
	   !grepl("^\\d+$", gene)) |>
    mutate(cell_type = factor(cell_type, levels = c("Naive B", "USM B", "SM B", "DN B", "Plasmablast")))

#pval_df <- 
#    full_join(dis_activ_b, dis_state_b, 
#	  join_by(cell_type, gene), 
#	  suffix = c("_activ", "_state")) |>
#    select(cell_type, gene, pvalue_activ, pvalue_state)
#
#
#pval_plot <- 
#    ggplot(pval_df, aes(x = pvalue_state, y = pvalue_activ)) +
#    geom_abline() +
#    geom_point(size = .2, alpha = .2) +
#    scale_x_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
#    scale_y_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
#    facet_wrap(~cell_type, nrow = 1) +
#    theme_bw()
#
#ggsave("./pval_plot.png", pval_plot, width = 8, height = 2)


perez_state <- 
    expand_grid(perez, distinct(dis_state_b_all, cell_type)) |>
    left_join(dis_state_b_all, join_by(gene_name == gene, cell_type)) |>
    mutate(category = case_when(log_fc > 0 & fdr <= 0.05 ~ "up",
				log_fc < 0 & fdr <= 0.05 ~ "down",
				fdr > 0.05 ~ "n.s.",
				is.na(fdr) ~ "N.A.",
				TRUE ~ NA_character_)) |>
    count(cell_type, category) |> 
    mutate(category = fct_inorder(category))

perez_activ <- 
    expand_grid(perez, distinct(dis_activ_b_all, cell_type)) |>
    left_join(dis_activ_b_all, join_by(gene_name == gene, cell_type)) |>
    mutate(category = case_when(log_fc > 0 & fdr <= 0.05 ~ "up",
				log_fc < 0 & fdr <= 0.05 ~ "down",
				fdr > 0.05 ~ "n.s.",
				is.na(fdr) ~ "N.A.",
				TRUE ~ NA_character_)) |>
    count(cell_type, category) |>
    mutate(category = fct_inorder(category))

perez_state_plot <- 
    ggplot(perez_state) +
    geom_col(aes(x = n, y = cell_type, fill = category), 
             color = "black", position = "fill", linewidth = .25) +
    geom_vline(xintercept = .5, linetype = 2, linewidth = 0.25, color = "black") +
    scale_x_continuous(breaks = c(0, .5, 1)) +
    scale_fill_manual(values = c("N.A." = "grey", "down" = "001260", "n.s." = "white", "up" = "#AA4613")) +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          strip.text.x = element_text(size = 7),
          strip.text.y = element_text(size = 7),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9),
	  plot.title = element_text(size = 9)
    ) +
    labs(x = "Proportion of gene signature", y = NULL,
	 fill = "DE direction:",
	 title = "Perez et al. B cell eGenes in Nakano et al disease-state data")

perez_activ_plot <- 
    ggplot(perez_activ) +
    geom_col(aes(x = n, y = cell_type, fill = category), 
             color = "black", position = "fill", linewidth = .25) +
    geom_vline(xintercept = .5, linetype = 2, linewidth = 0.25, color = "black") +
    scale_x_continuous(breaks = c(0, .5, 1)) +
    scale_fill_manual(values = c("N.A." = "grey", "down" = "001260", "n.s." = "white", "up" = "#AA4613")) +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          strip.text.x = element_text(size = 7),
          strip.text.y = element_text(size = 7),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9),
	  plot.title = element_text(size = 9)
    ) +
    labs(x = "Proportion of gene signature", y = NULL,
	 fill = "DE direction:",
	 title = "Perez et al. B cell eGenes in Nakano et al disease-activity data")

ggsave("perez_egenes_in_nakano.png",
       plot_grid(perez_state_plot, perez_activ_plot, ncol = 1),
       width = 6, height = 4)





dis_activ_b |> 
    filter(gene %in% tfs_curated$Gene, fdr < 0.05) |>
    arrange(cell_type, sign(log_fc), fdr) |>
    print(n = Inf)

dis_state_b |> 
    filter(gene %in% tfs_curated$Gene, fdr < 0.05) |>
    arrange(cell_type, sign(log_fc), fdr) |>
    print(n = Inf)


################################################################################
# Volcano plot for IRF5
library(ggrepel)

cell_colors <- 
    c("Naive CD4" = "#7E3F3B",
      "Mem CD4" = "#F0948F",
      "Th1" = "#B33C3A",
      "Th2" = "#EA4439",
      "Th17" = "#D6A08A",
      "Tfh" = "#D16C61",
      "Fr. I nTreg" = "#F2A767",
      "Fr. II eTreg" = "#E67744",
      "Fr. III T" = "#F6D9B6",
      "Naive CD8" = "#875D42",
      "CM CD8" = "#BB9062",
      "EM CD8" = "#E0B27C",
      "TEMRA CD8" = "#AB9082",
      "NK" = "#D6CB57",
      "Naive B" = "#39593F",
      "USM B" = "#528045",
      "SM B" = "#92B54F",
      "DN B" = "#72D15B",
      "Plasmablast" = "#B1E44D",
      "CL Mono" = "#5F619A",
      "CD16p Mono" = "#6DA3F9",
      "Int Mono" = "#B5C5F8",
      "NC Mono" = "#7C96F1", 
      "mDC" = "#72D7EE",
      "pDC" = "#BBF9F9",
      "Neu" = "#CCCDCD",
      "LDG" = "#A2A2A2"
    )

dis_state_irf5 <- 
    dis_state |> 
    filter(gene == "IRF5") |> 
    mutate(cell_type = factor(cell_type, levels = names(cell_colors)))

fdr_thres <- 
    dis_state |>
    filter(fdr <= 0.05) |>
    slice_max(pvalue, with_ties = FALSE) |>
    pull(pvalue)


volcano_plot <- 
    ggplot(data = dis_state_irf5) + 
    geom_hline(yintercept = -log10(fdr_thres), linetype = 2, linewidth = .2, color = "red") +
    geom_point(aes(x = log_fc, y = -log10(pvalue), fill = cell_type),
	       shape = 21, stroke = .2, size = 3) +
    geom_label_repel(data = filter(dis_state_irf5, log_fc >= 0.25, fdr <= 0.05),
		     aes(x = log_fc, y = -log10(pvalue), label = cell_type, color = cell_type),
		     size = 2, fill = alpha("white", .5),
		     box.padding = .3, 
		     point.padding = .25,
		     label.padding = unit(0.1, "lines"),
		     min.segment.length = 0, 
		     segment.size = .25, 
		     max.overlaps = Inf, 
		     segment.color = "black",
		     fontface = "italic") +
    scale_x_continuous(limits = c(-1.5, 1.5)) +
    scale_color_manual(values = cell_colors) +
    scale_fill_manual(values = cell_colors) +
    theme_minimal() +
    theme(
	  axis.title = element_text(size = 8),
	  axis.text = element_text(size = 7),
	  panel.grid.major.x = element_line(color = "gray90"),
	  panel.grid.major.y = element_line(color = "gray90"),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")
    ) + 
    labs(x = "log(FC)", y = expression(-log[10](p))) +
    guides(color = "none", fill = "none")

ggsave("./IRF5_volcano.png", volcano_plot, width = 3, height = 2, dpi = 300)

volcano_plot <- 
    ggplot(data = dis_state_irf5) + 
    geom_hline(yintercept = -log10(fdr_thres), linetype = 2, linewidth = .2, color = "red") +
    geom_point(aes(x = log_fc, y = -log10(pvalue), fill = cell_type),
	       shape = 21, stroke = .2, size = 3) +
    geom_label_repel(data = filter(dis_state_irf5, log_fc >= 0.25, fdr <= 0.05),
		     aes(x = log_fc, y = -log10(pvalue), label = cell_type, color = cell_type),
		     size = 2, fill = alpha("white", .5),
		     box.padding = .3, 
		     point.padding = .25,
		     label.padding = unit(0.1, "lines"),
		     min.segment.length = 0, 
		     segment.size = .25, 
		     max.overlaps = Inf, 
		     segment.color = "black",
		     fontface = "italic") +
    scale_x_continuous(limits = c(-0.5, 1.5)) +
    scale_color_manual(values = cell_colors) +
    scale_fill_manual(values = cell_colors) +
    theme_minimal() +
    theme(
	  axis.title = element_text(size = 8),
	  axis.text = element_text(size = 7),
	  legend.title = element_text(size = 8),
	  legend.text = element_text(size = 7, margin = margin(l = 0)),
	  legend.key.height = unit(1, "null"),
	  panel.grid.major.x = element_line(color = "gray90"),
	  panel.grid.major.y = element_line(color = "gray90"),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")
    ) + 
    labs(x = "log(FC)", y = expression(-log[10](p)), fill = "Cell type:") +
    guides(color = "none")

ggsave("./IRF5_volcano.png", volcano_plot, width = 5, height = 2, dpi = 300)



#library(JASPAR2024)
#library(TFBSTools)
#
#pfm_list <- getMatrixSet(db(JASPAR2024()), list(species = 9606))
#
#
#tf_families <- tags(pfm_list) |> map("family") |> map(~paste(., collapse = "::")) |> enframe("motif_id", "tf_family") |> unnest(cols = tf_family)
#tf_names <- map_chr(pfm_list, "name") |> enframe("motif_id", "tf_name")
#
#tf_aliases <- tags(pfm_list) |> map("alias") |> enframe("motif_id", "alias") |> unnest(cols = alias)
#
#pfm_df <- left_join(tf_names, tf_families)
#
#library(clusterProfiler)
#
#?GSEA
#
#gene_set_df <- 
#    gene_sets |> 
#    map(function(x) tibble(gene = x)) |>
#    bind_rows(.id = "TERM")
#
#gsea_res <- 
#    GSEA(
#	 geneList = sort(bcell_gene_lists[[2]], decreasing = TRUE),
#	 TERM2GENE = gene_set_df,
#	 pvalueCutoff = 1,
#	 nPerm = 1000,
#	 by = "DOSE",
#	 verbose = TRUE
#    )
#
#as.data.frame(gsea_res) |> as_tibble()
#
#
#



#library(piano)
#
#gsc <- 
#    gene_sets |> 
#    map(function(x) tibble(g = x)) |>
#    bind_rows(.id = "s") |>
#    select(g, s) |>
#    loadGSC(type = "data.frame")
#
#runGSA(
#       geneLevelStats = bcell_gene_lists[[11]],
#       gsc            = gsc,
#       geneSetStat    = "gsea",          
#       signifMethod   = "geneSampling",  
#       nPerm          = 10000, 
#       adjMethod      = "fdr",
#       ncpus          = 1,               
#       verbose        = TRUE)
#
#daemons(16)
#
#gsa_res <- 
#    map(bcell_gene_lists[11],
#	in_parallel(function(x) piano::runGSA(
#				  geneLevelStats = x,
#				  gsc            = gsc,
#				  geneSetStat    = "gsea",          
#				  signifMethod   = "geneSampling",  
#				  nPerm          = 10000, 
#				  adjMethod      = "fdr",
#				  ncpus          = 1,               
#				  verbose        = TRUE),
#		    gsc = gsc))
#
#
#daemons(0)
#
#gsa_res_df <- 
#    map_dfr(gsa_res, GSAsummaryTable, .id = "contrast") |>
#    janitor::clean_names()
#
