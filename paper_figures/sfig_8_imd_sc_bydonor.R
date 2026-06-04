library(tidyverse)
library(Seurat)

features_df <-
    "../04_citeseq/data/cellranger/1984/filtered_feature_bc_matrix/features.tsv.gz" |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

bcells <- read_rds("../04_citeseq/2-processing/data/v4_seurat_qced.rds")

marker_genes <- c("CTSH", "CXCR5", "IRF5", "IRF8", "MYC", "RGS1", "ZBTB38")

marker_genes_df <-
    features_df |>
    filter(phenotype == "Gene Expression",
           gene_name %in% marker_genes) |>
    select(gene_id, gene_name) |>
    mutate(gene_name = factor(gene_name, levels = marker_genes)) |>
    arrange(gene_name)

markers_plot_data <- 
    DotPlot(object = bcells, 
	    features = marker_genes_df$gene_id, 
	    split.by = 'donor_id',
	    cols = rep("black", 5)) |>
    {function(x) as_tibble(x$data)}() |>
    left_join(marker_genes_df, join_by(features.plot == gene_id)) |>
    separate(id, c("cluster", "donor"), sep = "_") |>
    select(cluster, donor, gene_name, avg.exp.scaled, pct.exp) |>
    mutate(cluster = factor(cluster, levels = sort(as.numeric(unique(cluster)))),
           donor = factor(donor, labels = paste("Donor", 1:5))) |>
    arrange(donor, gene_name, cluster)

dotplot_markers <- 
    ggplot(markers_plot_data, aes(x = donor, y = gene_name)) +
    geom_point(aes(size = pct.exp, fill = avg.exp.scaled), stroke = 0.2, shape = 21) +
    scale_fill_gradient2(low = "Light Sky Blue", mid = "lightyellow", high = "Dark Red",
                         midpoint = 0) +
    scale_size(range = c(0.1, 2.5), 
               limits = c(0, 100),
               breaks = c(0, 33, 66, 100)) +
    facet_wrap(~cluster, nrow = 1) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, face = 'italic'),
        panel.grid = element_line(linewidth = .25, color = "grey90"),
        panel.spacing.x = unit(.1, "lines"),
        legend.margin = margin(t = 0, r = 0.2, b = 0, l = 0, unit = "lines"),
        legend.title = element_text(size = 8),
        legend.key.spacing.y = unit(0, "lines"),
        legend.box.margin = margin(l = 0, b = 5),
        strip.clip = "off",
        plot.margin = margin(0, 0.1, 0, 0, unit = "lines"),
        plot.background = element_rect(fill = "white", color = "white")
    ) +
    guides(fill = guide_colorbar(order = 1, position = "right",
                                 barwidth = .25, barheight = 3)) +
    labs(x = NULL, y = NULL, fill = "Scaled\nExpression", size = "%\nExpressed")

ggsave("./sfigs/sfig8_imd_sc_bydonor.png", dotplot_markers, width = 6.5, height = 2)















