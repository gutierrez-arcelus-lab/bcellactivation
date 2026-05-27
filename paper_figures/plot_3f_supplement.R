library(tidyverse)

cellranger_dir_1984 <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq",
              "SN0268787/broad/hptmp/curtism/bwh10x/KW10456_Maria",
              "221101_10X_KW10456_bcl/cellranger-6.1.1/GRCh38/BRI-1984/outs",
              "filtered_feature_bc_matrix")

features_df <- 
    file.path(cellranger_dir_1984, "features.tsv.gz") |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

bcells <- read_rds("../citeseq/data/seuratv4_qced.rds")

marker_genes <- c("CTSH", "CXCR5", "IRF5", "IRF8", "MYC", "RGS1", "ZBTB38")

marker_genes_df <-
    features_df |>
    filter(phenotype == "Gene Expression",
           gene_name %in% marker_genes) |>
    select(gene_id, gene_name) |>
    mutate(gene_name = factor(gene_name, levels = marker_genes)) |>
    arrange(gene_name)


markers_plot_data <- 
    Seurat::DotPlot(object = bcells, 
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

ggsave("./fig3f_supplement.png", dotplot_markers, width = 6.5, height = 2)















