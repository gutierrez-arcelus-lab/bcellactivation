library(tidyverse)
library(readxl)
#library(edgeR)
library(fgsea)
library(mirai)

# In Nakano et al, they say they followed analysis from Ota et al 2021, where they used Gencode v27
# Still, some gene names don't match those annotated in Gencode v27
gene_ids <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
              "gencode.v27.primary_assembly.annotation.gtf.gz") |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X3 == "gene") |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
           gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(gene_id, gene_name) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

dis_state <- 
    read_excel("./data/mmc2.xlsx", sheet = 2, skip = 2) |>
    janitor::clean_names()

gene_sets <- 
    dis_state |> 
    filter(cell_type %in% c("Naive B", "USM B", "SM B", "DN B", "Plasmablast"),
	   log_fc > 0, fdr < 0.01) |>
    left_join(gene_ids, join_by(gene == gene_name)) |>
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
library(glue)
library(cowplot)

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

gene_sets_df <- 
    enframe(gene_sets, "cell_type", "gene_id") |>
    unnest(cols = c(gene_id))

vs_unstim0_plot_data <- 
    gene_sets_df |>
    expand_grid(distinct(vs_unstim0, condition, stim1, t1)) |>
    left_join(vs_unstim0, join_by(gene_id, condition, stim1, t1)) |>
    mutate(category = case_when(logFC > 0 & FDR < 0.05 ~ "up",
				logFC < 0 & FDR < 0.05 ~ "down",
				FDR > 0.05 ~ "n.s.",
				is.na(FDR) ~ "N.A.",
				TRUE ~ NA_character_)) |>
    count(cell_type, condition, stim1, t1, category) |>
    mutate(category = fct_inorder(category),
           cell_type = factor(cell_type, levels = c("Naive B", "USM B", "SM B", "DN B", "Plasmablast")))


vs_unstim0_plot <- 
    ggplot(vs_unstim0_plot_data) +
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
    labs(x = "Hours", y = "Proportion of gene signature", fill = "DE direction:")

ggsave("./vs_unstim0.png", vs_unstim0_plot, width = 8, height = 4)







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
