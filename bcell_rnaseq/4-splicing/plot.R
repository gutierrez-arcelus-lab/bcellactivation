library(tidyverse)
library(glue)
library(ggrepel)

read_leaf <- function(stim) {
    
    sig <- 
        glue("./results/{stim}_cluster_significance.txt") |>
        read_tsv() |>
	filter(status == "Success") |>
        separate(cluster, c("chr", "cluster"), sep = ":") |>
        select(cluster, p, padj = p.adjust, genes)

    eff <- 
        glue("./results/{stim}_effect_sizes.txt") |>
        read_tsv() |>
        mutate(cluster = sub("^.*(clu.*)$", "\\1", intron)) |>
        select(cluster, logef, deltapsi)

    inner_join(sig, eff, join_by(cluster))
}

leaf_df <- 
    c("TLR7", "BCR", "DN2") |> 
    {function(x) setNames(x, x)}() |>
    map_dfr(~read_leaf(.x) |>
	    group_by(cluster) |>
	    slice_min(p) |>
	    slice_max(abs(deltapsi)) |>
	    ungroup(),
	    .id = "stim") |>
    mutate(stim = fct_inorder(stim))

plot_labels <- 
    leaf_df |>
    filter(padj <= 0.05) |>
    group_by(stim) |>
    top_n(20, -log10(p)) |>
    ungroup() |>
    mutate(genes = str_trunc(genes, 10))

leaf_plot <- 
    ggplot(leaf_df, aes(x = deltapsi, y = -log10(p))) +
    geom_point(size = .5) +
    geom_text_repel(data = plot_labels, 
		    aes(x = deltapsi, y = -log10(p), label = genes),
		    size = 2,
		    min.segment.length = 0,
		    segment.size = .25,
		    max.overlaps = Inf) +
    facet_wrap(~stim, ncol = 1) +
    theme_bw() +
    theme(panel.grid = element_blank())

ggsave("./plots/leafcutter_ds.png", leaf_plot, height = 8, width = 5)
