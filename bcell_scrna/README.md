CITE-seq Pilots
================

Raw counts
----------

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

QC
--

Here we use the miQC package to model the percentage of mitochondrial
reads and number of genes, in order to identify and remove compromised
cells.

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

    # Removing 4104 out of 13946 cells.

    # Removing 723 out of 10562 cells.

Demultiplex cells based on HTO
------------------------------

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

In the plot above, we see 2 differences between pilots 1 and 2:

-   In pilot 1, HTO counts in higher, including for those cells being
    classified as “Negatives”;

-   In pilot 2, HTO counts for doublets not that higher than those for
    singlets.

Regarding the last point, doublets don’t necessarily have higher counts,
but a mixture of HTOs.

However, a lot of droplets being classified as doublets in pilot 2 have
a profile resembling singlets in respect to the mixture of HTO, as we
can see in the plots below.

    plot_admix <- function(seurat_obj, stim_colors) {

        meta_df <- seurat_obj@meta.data %>%
        as_tibble(rownames = "barcode") %>%
        select(barcode, stim = HTO_maxID, hto_class = HTO_classification.global)
      
        hto_df <- seurat_obj@assays$HTO@counts %>%
        as_tibble(rownames = "hto") %>%
        pivot_longer(-hto, names_to = "barcode") %>%
        left_join(meta_df)

        hto_top <- hto_df %>%
        group_by(barcode, hto_class) %>%
        slice_max(n = 1, order_by = value) %>%
        select(barcode, hto_class, top_hto = hto) %>%
        ungroup()

        hto_plot_df <- hto_df %>%
        left_join(hto_top, by = c("barcode", "hto_class")) %>%
        mutate_at(vars(top_hto, stim), ~factor(., levels = names(stim_colors)))

        ggplot(hto_plot_df, aes(reorder_within(barcode, by = value, within = stim), value)) +
        geom_col(aes(fill = hto), position = "fill", width = 1.01, show.legend = FALSE) +
        scale_fill_manual(values = stim_colors) +
        facet_grid(~stim, scales = "free_x", space = "free") +
        theme_bw() +
        theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          strip.text = element_blank()) +
        labs(x = NULL, y = NULL)
    }

    Idents(bcells_1_m1_qc_demux) <- "HTO_classification.global"
    Idents(bcells_1_m2_qc_demux) <- "HTO_classification.global"
    Idents(bcells_2_m1_qc_demux) <- "HTO_classification.global"
    Idents(bcells_2_m2_qc_demux) <- "HTO_classification.global"


    admix_1 <- plot_admix(subset(bcells_1_m1_qc_demux, idents = "Singlet"), pilot1_pal) +
        labs(title = "Pilot 1: Singlets", subtitle = "margin = 1")

    admix_2 <- plot_admix(subset(bcells_1_m2_qc_demux, idents = "Singlet"), pilot1_pal) +
        labs(title = "Pilot 1: Singlets", subtitle = "margin = 2")

    admix_3 <- plot_admix(subset(bcells_1_m1_qc_demux, idents = "Doublet"), pilot1_pal) +
        labs(title = "Pilot 1: Doublets", subtitle = "margin = 1")

    admix_4 <- plot_admix(subset(bcells_1_m2_qc_demux, idents = "Doublet"), pilot1_pal) +
        labs(title = "Pilot 1: Doublets", subtitle = "margin = 2")

    admix_5 <- plot_admix(subset(bcells_2_m1_qc_demux, idents = "Singlet"), pilot2_pal) +
        labs(title = "Pilot 2: Singlets", subtitle = "margin = 1")

    admix_6 <- plot_admix(subset(bcells_2_m2_qc_demux, idents = "Singlet"), pilot2_pal) +
        labs(title = "Pilot 2: Singlets", subtitle = "margin = 2")

    admix_7 <- plot_admix(subset(bcells_2_m1_qc_demux, idents = "Doublet"), pilot2_pal) +
        labs(title = "Pilot 2: Doublets", subtitle = "margin = 1")

    admix_8 <- plot_admix(subset(bcells_2_m2_qc_demux, idents = "Doublet"), pilot2_pal) +
        labs(title = "Pilot 2: Doublets", subtitle = "margin = 2")

    plot_grid(admix_1, admix_2, admix_3, admix_4, admix_5, admix_6, admix_7, admix_8, ncol = 1)

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

So, Seurat might be doing a good job in classifying Singlets. For
exemple, see the plot below for Pilot 1, margin = 2.

    hto_counts <- bcells_1_m2_qc_demux@assays$HTO@counts %>%
        as_tibble(rownames = "hto") %>%
        pivot_longer(-hto, names_to = "barcode", values_to = "count") %>%
        left_join(as_tibble(bcells_1_m2_qc_demux@meta.data, rownames = "barcode") %>%
                  select(barcode, hto_max = HTO_maxID, hto_class = HTO_classification.global)) %>%
        mutate_at(vars(hto_max, hto), ~factor(., levels = pilot1_stim)) %>%
        mutate(hto_class = factor(hto_class, levels = c("Singlet", "Doublet", "Negative"))) %>%
        arrange(hto_class, hto, hto_max) %>%
        mutate(logcount = ifelse(count > 0, log10(count), 0))

    hto_counts %>%
        filter(hto_class == "Singlet") %>%
        ggplot(aes(x = hto, y = logcount)) +
        geom_violin(aes(fill = hto), size = .25) +
        geom_quasirandom(data = filter(hto_counts, hto_class == "Negative"),
                         method = "smiley", groupOnX = TRUE,
                         shape = 21, fill = NA, size = .5, width = .25) +
        scale_fill_manual(values = pilot1_pal) +
        facet_wrap(~hto_max, nrow = 2) +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank()) +
        labs(y = "Log10(count)",
             fill = "HTO", 
             title = "Distribution of all HTO counts for singlets in each condition",
             subtitle = "Points represent droplets classified as Negatives")

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Let’s take a look at the clustering performed by Seurat to classify
cells.

    plot_clusters <- function(seurat_obj, stim_colors) {
        
        DefaultAssay(seurat_obj) <- "HTO"
        
        seurat_obj <- seurat_obj %>%
        ScaleData(., features = rownames(.), verbose = FALSE) %>%
        RunPCA(., features = rownames(.), approx = FALSE, verbose = FALSE)

        hto_k <- clara(t(seurat_obj@assays$HTO@data), nrow(seurat_obj) + 1, samples = 100)
        kmeans_df <- tibble(barcode = names(hto_k$clustering),
                cluster = hto_k$clustering)

        p1 <- seurat_obj@reductions$pca@cell.embeddings %>%
        as_tibble(rownames = "barcode") %>%
        left_join(kmeans_df) %>%
        ggplot(aes(PC_1, PC_2, color = factor(cluster))) +
        geom_point(size = .75) +
        scale_color_manual(values = pals::kelly(n = nrow(seurat_obj) + 1)) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        guides(color = guide_legend(override.aes = list(size = 2))) +
        labs(color = "Cluster")

        p2 <- seurat_obj@reductions$pca@cell.embeddings %>%
        as_tibble(rownames = "barcode") %>%
        left_join(as_tibble(seurat_obj@meta.data, rownames = "barcode") %>%
              select(barcode, hto_class = HTO_classification) %>%
              mutate(hto_class = ifelse(grepl("_", hto_class), "Doublet", hto_class))) %>%
        ggplot(aes(PC_1, PC_2, color = hto_class)) +
        geom_point(size = .75) +
        scale_color_manual(values = c(stim_colors, "Negative" = "red")) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        guides(color = guide_legend(override.aes = list(size = 2))) +
        labs(color = "HTO\nclassification")

        plot_grid(p1, p2, nrow = 1)
    }

    plot_grid(plot_clusters(bcells_1_m1_qc_demux, pilot1_pal) +
          labs(title = "Pilot 1, margin = 1"),
          plot_clusters(bcells_1_m2_qc_demux, pilot1_pal) +
          labs(title = "Pilot 1, margin = 2"),
          plot_clusters(bcells_2_m1_qc_demux, pilot2_pal) +
          labs(title = "Pilot 2, margin = 1"),
          plot_clusters(bcells_2_m2_qc_demux, pilot2_pal) +
          labs(title = "Pilot 2, margin = 2"),
          ncol = 1)

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Alternative demultiplex algorithm
---------------------------------

-   Cells with total HTO counts &lt; 20 times the number of different
    hashtags = “Negatives”;

-   Cells with a dominant hashtag (75% of total counts) = “Singlets”;

-   Cells with a mixture = “Doublets”

<!-- -->

    custom_class_pilot1 <- bcells_1_m1_qc_demux@assays$HTO@counts %>%
        as_tibble(rownames = "hto") %>%
        pivot_longer(-hto, names_to = "barcode", values_to = "count") %>%
        group_by(barcode) %>%
        mutate(pct = count/sum(count)) %>%
        summarise(hto_class = case_when(sum(count) < n_distinct(hto) * 20 ~ "Negative",
                        sum(count) >= n_distinct(hto) * 20 & any(pct >= .75) ~ "Singlet",
                        sum(count) >= n_distinct(hto) * 20 & !any(pct >= .75) ~ "Doublet",
                        TRUE ~ NA_character_)) %>%
        ungroup()

    custom_class_pilot2 <- bcells_2_m1_qc_demux@assays$HTO@counts %>%
        as_tibble(rownames = "hto") %>%
        pivot_longer(-hto, names_to = "barcode", values_to = "count") %>%
        group_by(barcode) %>%
        mutate(pct = count/sum(count)) %>%
        summarise(hto_class = case_when(sum(count) < n_distinct(hto) * 20 ~ "Negative",
                        sum(count) >= n_distinct(hto) * 20 & any(pct >= .75) ~ "Singlet",
                        sum(count) >= n_distinct(hto) * 20 & !any(pct >= .75) ~ "Doublet",
                        TRUE ~ NA_character_)) %>%
        ungroup()

    write_tsv(custom_class_pilot1, "./custom_class_pilot1.tsv")
    write_tsv(custom_class_pilot2, "./custom_class_pilot2.tsv")

Comparison of the number of cells after QC:

    list(
         "Pilot 1, margin=1" = count(meta_1_1, hto_class),
         "Pilot 1, margin=2" = count(meta_1_2, hto_class),
         "Pilot 2, margin=1" = count(meta_2_1, hto_class),
         "Pilot 2, margin=2" = count(meta_2_2, hto_class),
         "Pilot 1, Alternative" = count(custom_class_pilot1, hto_class),
         "Pilot 2, Alternative" = count(custom_class_pilot2, hto_class)) %>%
    bind_rows(.id = "set") %>%
    pivot_wider(names_from = hto_class, values_from = n)

    # # A tibble: 6 × 4
    #   set                  Doublet Negative Singlet
    #   <chr>                  <int>    <int>   <int>
    # 1 Pilot 1, margin=1       3023     1115    5654
    # 2 Pilot 1, margin=2       3143      131    6518
    # 3 Pilot 2, margin=1       3961       40    5814
    # 4 Pilot 2, margin=2       3914       33    5868
    # 5 Pilot 1, Alternative    2825        5    6962
    # 6 Pilot 2, Alternative     985      407    8423
