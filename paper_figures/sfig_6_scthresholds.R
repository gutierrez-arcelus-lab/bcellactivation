library(tidyverse)
library(Seurat)
library(ggridges)
library(patchwork)

stim_colors <- 
    read_tsv("./figure_colors_final.txt") |>
    unite("stim", c(Condition, Time), sep = " ") |>
    mutate(stim = paste0(stim, "h")) |>
    select(stim, Hex) |>
    deframe()

# data 
bcells <- 
    read_rds("../04_citeseq/2-processing/data/v4_seurat_qced.rds")


bcells@meta.data <- 
    bcells@meta.data |>
    mutate(dmm_hto_call = str_replace(dmm_hto_call, "IL4", "IL-4c"),
	   dmm_hto_call = str_replace(dmm_hto_call, "TLR7", "TLR7c"),
	   dmm_hto_call = str_replace(dmm_hto_call, "BCR", "BCRc"),
	   dmm_hto_call = str_replace(dmm_hto_call, "DN2", "DN2c"))

cluster_markers <- 
    read_tsv("../04_citeseq/2-processing/data/v4_cluster_markers.tsv")

adt_df <- 
    bcells@assays$ADT@data[c('IgD', 'CD27', 'CD11c'), ] |>
    t() |>
    as_tibble(rownames = 'barcode') |>
    left_join(bcells@meta.data |> as_tibble(rownames = 'barcode') |> select(barcode, dmm_hto_call, seurat_clusters)) |>
    select(barcode, hto = dmm_hto_call, cluster = seurat_clusters, everything())

gates_genes <- 
    cluster_markers |> 
    filter(gene_name %in% c('MKI67', 'XBP1')) |>
    distinct(gene, gene_name)

rna_df <-
    bcells@assays$RNA@data[gates_genes$gene, ] |>
    as.matrix() |>
    t() |>
    as_tibble(rownames = 'barcode') |>
    pivot_longer(-barcode, names_to = 'gene') |>
    left_join(gates_genes) |>
    select(-gene) |>
    pivot_wider(names_from = gene_name, values_from = value) |>
    left_join(bcells@meta.data |> as_tibble(rownames = 'barcode') |> select(barcode, dmm_hto_call, seurat_clusters)) |>
    select(barcode, hto = dmm_hto_call, cluster = seurat_clusters, everything())

threshold_df <-
   tribble(~marker, ~thr,
    "IgD", .75,
    "CD27", .25,
    "CD11c", .75,
    "MKI67", 1,
    "XBP1", 2.5,
   ) |>
    mutate(marker = fct_inorder(marker))

thres_plot <-
    left_join(adt_df, rna_df, join_by(barcode, hto, cluster)) |>
    pivot_longer(-c(barcode, hto, cluster), names_to = "marker") |>
    mutate(marker = fct_inorder(marker)) |>
        ggplot(aes(y = cluster, x = value)) +
        geom_density_ridges(size = .2) +
        geom_vline(data = threshold_df, aes(xintercept = thr), linewidth = .2, color = "red") +
        facet_wrap(~marker, nrow = 1, scales = "free_x") +
        theme_bw() +
        theme(panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_line(color = "grey96")
        )

igd_cd27_plot_density <-
   ggplot(adt_df, aes(x = IgD, y = CD27)) +
   scale_fill_continuous(type = "viridis") +
   geom_vline(xintercept = .75, color = "red", linewidth = .2) +
   geom_hline(yintercept = .25, color = "red", linewidth = .2) +
   facet_wrap(~factor(hto, levels = names(stim_colors)), nrow = 1) +
   geom_hex(aes(fill = after_stat(count/max(count))), bins = 100) +
   theme_bw() +
   theme(panel.grid = element_blank(),
         strip.text = element_text(size = 8),
         panel.spacing.x = unit(.1, "lines")
         ) +
   guides(fill = "none")

ggsave("./sfigs/sfig6_scthreshold.png",
      thres_plot / igd_cd27_plot_density + plot_layout(heights = c(1, 0.45)) + plot_annotation(tag_levels = "A"),
      width = 6.5, height = 4)
