library(tidyverse)
library(glue)
library(readxl)

################################################################################
# Low-input data
################################################################################
stim_colors <- 
    read_tsv("./figure_colors_final.txt") |>
    unite("stim", c(Condition, Time), sep = "_") |>
    mutate(stim = paste0(stim, "hrs")) |>
    select(stim, Hex) |>
    deframe()

all_stims <-
    names(stim_colors) |>
    str_split("_") |>
    map_chr(1) |>
    unique()

diff_expr <- 
    "../01_rnaseq_lowinput/2_dge/results/edger/diff_expr_all_times_all_genes.tsv" |>
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
# Gencode annotations
################################################################################

# In Nakano et al, they say they followed analysis from Ota et al 2021, where they used Gencode v27
# Still, some gene names don't match those annotated in Gencode v27
gencode_v27 <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
              "gencode.v27.primary_assembly.annotation.gtf.gz") |>
    rtracklayer::import(feature.type = "gene") |>
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

ggsave("./sfigs/sfig5_nakano_sigs.png", nakano_plot, width = 6.5, height = 3.5)


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

