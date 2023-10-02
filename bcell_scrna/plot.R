library(Seurat)
library(tidyverse)
library(ggrepel)

genes_df <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/CITEseq_pilot/SN0231064/KW9100_Maria/210726_10X_KW9100-2_bcl/cellranger-6.0.1/GRCh38/BRI-1283/outs/filtered_feature_bc_matrix/features.tsv.gz" %>%
    read_tsv(col_names = c("gene", "gene_name", "phenotype")) %>%
    filter(phenotype == "Gene Expression") %>%
    select(gene, gene_name)

pilot1 <- readRDS("data/pilot1.Rdata")
pilot2 <- readRDS("data/pilot2.Rdata")

Idents(pilot1) <- "HTO_maxID"
Idents(pilot2) <- "HTO_maxID"

# pilot 1
res24_1 <- 
    FindMarkers(pilot1, 
                   ident.1 = "Res 24hr",
                   ident.2 = "Res 0hr",
                   only.pos = TRUE,
                   min.pct = 0.1,
                   logfc.threshold = 0) %>%
    rownames_to_column("gene") %>%
    as_tibble()

bcr24_1 <- 
    FindMarkers(pilot1, 
                   ident.1 = "BCR 24hr",
                   ident.2 = "Res 0hr",
                   only.pos = TRUE,
                   min.pct = 0.1,
                   logfc.threshold = 0) %>%
    rownames_to_column("gene") %>%
    as_tibble()

tlr24_1 <- 
    FindMarkers(pilot1, 
		ident.1 = "TLR7 24hr",
		ident.2 = "Res 0hr",
		only.pos = TRUE,
		min.pct = 0.1,
		logfc.threshold = 0) %>%
    rownames_to_column("gene") %>%
    as_tibble()


markers_1 <-
  bind_rows("Null" = res24_1,
	    "BCR" = bcr24_1,
            "TLR7" = tlr24_1,
            .id = "stim")

# pilot 2

il4_2 <- 
    FindMarkers(pilot2, 
		ident.1 = "IL4 24hr",
		ident.2 = "day0",
		only.pos = TRUE,
		min.pct = 0.1,
		logfc.threshold = 0) %>%
    rownames_to_column("gene") %>%
    as_tibble()

bcr24_2 <- 
    FindMarkers(pilot2, 
		ident.1 = "BCR 24hr",
		ident.2 = "day0",
		only.pos = TRUE,
		min.pct = 0.1,
		logfc.threshold = 0) %>%
    rownames_to_column("gene") %>%
    as_tibble()

tlr24_2 <- 
    FindMarkers(pilot2, 
		ident.1 = "TLR7 24hr",
		ident.2 = "day0",
		only.pos = TRUE,
		min.pct = 0.1,
		logfc.threshold = 0) %>%
    rownames_to_column("gene") %>%
    as_tibble()

markers_2 <-
  bind_rows("Null" = il4_2,
	    "BCR" = bcr24_2,
            "TLR7" = tlr24_2,
            .id = "stim")

# plot
plot_df <- 
    full_join(select(markers_1, stim, gene, eff_1 = avg_log2FC),
	      select(markers_2, stim, gene, eff_2 = avg_log2FC),
	      by = c("stim", "gene")) %>%
    mutate_at(vars(eff_1:eff_2), ~replace_na(., 0)) %>%
    left_join(genes_df, by = "gene") %>%
    select(stim, gene, gene_name, eff_1, eff_2) %>%
    mutate(stim = recode(stim, "Null" = "Unstim/IL4"),
	   stim = factor(stim, levels = c("Unstim/IL4", "BCR", "TLR7")))

out <- ggplot(plot_df, aes(eff_1, eff_2)) +
    geom_abline(linetype = 2) +
    geom_point(alpha = .25) +
    ggrepel::geom_text_repel(data = filter(plot_df, abs(eff_1 - eff_2) > 1.66),
			     aes(label = gene_name),
			     size = 1.5, segment.color = "grey") +
    facet_wrap(~stim, nrow = 1) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "Pilot 1", y = "Pilot 2",
	 title = "Log2(FC) of each stimulation at 24hr in respect to Unstim at day 0.")

ggsave("./plot.png", out, width = 6.5, height = 3)

# Other comparisons

bcr_1 <- FindMarkers(pilot2, 
		ident.1 = "BCR 24hr",
		ident.2 = "day0",
		only.pos = TRUE,
		min.pct = 0.1,
		logfc.threshold = 0) %>%
    rownames_to_column("gene") %>%
    as_tibble()

bcr_2 <- FindMarkers(pilot2, 
		ident.1 = "BCR 24hr",
		ident.2 = "IL4 24hr",
		only.pos = TRUE,
		min.pct = 0.1,
		logfc.threshold = 0) %>%
    rownames_to_column("gene") %>%
    as_tibble()

tlr7_1 <- FindMarkers(pilot2, 
		ident.1 = "TLR7 24hr",
		ident.2 = "day0",
		only.pos = TRUE,
		min.pct = 0.1,
		logfc.threshold = 0) %>%
    rownames_to_column("gene") %>%
    as_tibble()

tlr7_2 <- FindMarkers(pilot2, 
		ident.1 = "TLR7 24hr",
		ident.2 = "IL4 24hr",
		only.pos = TRUE,
		min.pct = 0.1,
		logfc.threshold = 0) %>%
    rownames_to_column("gene") %>%
    as_tibble()

plot_df2 <- 
    bind_rows("BCR + IL4 vs. day0" = bcr_1,
	      "BCR + IL4 vs. IL4 alone" = bcr_2,
	      "TLR7 + IL4 vs. day0" = tlr7_1,
	      "TLR7 + IL4 vs. IL4 alone" = tlr7_2,
	      .id = "stim") %>%
    left_join(genes_df, by = "gene") %>%
    select(stim, gene, gene_name, eff = avg_log2FC) %>%
    extract(stim, c("stim", "comparison"), "(BCR|TLR7) (\\+.+)") %>%
    pivot_wider(names_from = comparison, values_from = eff) %>%
    mutate_at(vars(4:5), ~replace_na(., 0))

out2 <- ggplot(plot_df2, aes(`+ IL4 vs. day0`, `+ IL4 vs. IL4 alone`)) +
    geom_abline(linetype = 2) +
    geom_point(alpha = .1) +
    geom_text(data = filter(plot_df2, abs(`+ IL4 vs. day0` - `+ IL4 vs. IL4 alone`) > 1),
	       aes(label = gene_name),
	       size = 1, fontface = "bold", check_overlap = TRUE) +
    facet_wrap(~stim) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(title = "Log2(FC) to evaluate the effect of IL4 in pilot 2.")

ggsave("./plot-2.png", out2, width = 6, height = 3)


