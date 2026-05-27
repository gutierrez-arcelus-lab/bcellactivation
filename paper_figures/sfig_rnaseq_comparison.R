library(tidyverse)
library(ggpointdensity)

dge_low <- 
    "../bcell_lowinput/results/edger/diff_expr_all_times_all_genes.tsv" |>
    read_tsv() |>
    filter(contrast %in% c("TLR7.24-Unstim.0", "BCR.24-Unstim.0", "DN2.24-Unstim.0")) |>
    mutate(stim = str_extract(contrast, "(.+)-.+", group = TRUE),
	   stim = str_remove(stim, "\\.\\d+"),
	   gene_id = str_remove(gene_id, "\\.\\d+")) |>
    select(stim, gene_id, logFC:FDR)

dge_high <- 
    "../bcell_rnaseq/5-dge/results_v38.tsv" |>
    read_tsv()

dge_merge <- 
    inner_join(dge_low |> filter(FDR <= 0.05) |> select(stim, gene_id, logFC),
	       select(dge_high, stim, gene_id, logFC),
	       join_by(stim, gene_id),
	       suffix = c("_low", "_high")) |>
    mutate(stim = recode(stim, "TLR7" = "TLR7c", "BCR" = "BCRc", "DN2" = "DN2c"),
	   stim = factor(stim, levels = c("TLR7c", "BCRc", "DN2c")))

cor_df <- 
    dge_merge |>
    group_by(stim) |>
    summarize(r = cor(logFC_high, logFC_low, method = "pearson")) |>
    ungroup() |>
    mutate(label = paste0("R = ", round(r, 2)))

dge_plot <-
    ggplot(dge_merge) +
    geom_pointdensity(aes(x = logFC_low, y = logFC_high),
		      size = .05, show.legend = FALSE) +
    scale_color_viridis_c() +
    geom_text(data = cor_df,
	      aes(x = -Inf, y = Inf, label = label),
	      hjust = -.2, vjust = 1.5, 
	      inherit.aes = FALSE,
	      size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~stim) +
    theme_bw() +
    theme(
	  axis.title = element_text(size = 9),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  strip.background = element_rect(fill = "grey90"),

    ) +
    labs(
	 x = "Low-input LFC (N=6)",
	 y = "Standard-input LFC (N=16)"
    )

ggsave("./sfig_rnaseq_comparison.png", width = 5, height = 2, dpi = 300)


select(dge_low, stim, gene_id) |> mutate(low = 1) |> count(stim)
select(dge_high, stim, gene_id) |> mutate(high = 1) |> count(stim)


gene_annot <-
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf.gz" |>
    rtracklayer::import(feature.type = "gene") |>
    as.data.frame() |>
    as_tibble() |>
    select(gene_id, gene_name, gene_type) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

dge_merge |> 
    filter(logFC_low < 0, logFC_high > 1) |>
    left_join(gene_annot) |>
    arrange(desc(abs(logFC_low - logFC_high))) |> count(stim)


    


