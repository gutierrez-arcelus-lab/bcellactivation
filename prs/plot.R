library(tidyverse)
library(ggbeeswarm)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(uwot)

sle_data <- 
    "./mgb_data/sle_data.tsv" |>
    read_tsv(col_types = c(.default = "c"))

batches <- sprintf("04%02d", c(1:8, 10))

# chrX homozygosity
chrx <- 
    sprintf("/temp_work/ch229163/vcf/prs/chrx/chrX.MGB.%s.het", batches) |>
    setNames(batches) |>
    map_dfr(~read_tsv(.) |>
	    extract(INDV, c("subject_id"), "[^-]+-(\\d+)", remove = FALSE) |>
	    mutate(hom = `O(HOM)`/N_SITES) |>
	    select(subject_id, hom) |>
	    inner_join(sle_data, by = "subject_id"),
	    .id = "batch")

chrx |> 
    filter(gender == "Female") |> 
    arrange(desc(hom))

chrx_plot_batches <- 
    ggplot(chrx, aes(x = batch, y = hom)) +
    geom_quasirandom(method = "smiley", alpha = .25, size = .25) +
    theme_bw() +
    theme(plot.background = element_rect(fill = "white", color = "white"),
	  panel.grid = element_blank(),
	  legend.position = "none",
	  plot.title = element_text(size = 10)) +
    labs(x = NULL, y = NULL, title = "Homozygosity at chrX per batch")

chrx_plot <- 
    ggplot(chrx, aes(x = gender, y = hom, color = gender)) +
    geom_quasirandom(method = "smiley", alpha = .25, size = .5) +
    scale_color_manual(values = c("Male" = "midnightblue", "Female" = "tomato4")) +
    facet_wrap(~group) +
    theme_bw() +
    theme(plot.background = element_rect(fill = "white", color = "white"),
	  panel.grid = element_blank(),
	  legend.position = "none",
	  plot.title = element_text(size = 10) ) +
    labs(x = NULL, y = NULL, title = 'Homozygosity at chrX per "gender"')

ggsave("./plots/chrx.png", chrx_plot_batches / chrx_plot, height = 4, width = 5)




## PCA on genotype data
### Metadata for 1000G samples
index_1000G <- 
    "https://raw.githubusercontent.com/igsr/1000Genomes_data_indexes/master/data_collections/1000_genomes_project/1000genomes.sequence.index"

sample_annotation <- 
    read_tsv(index_1000G, comment = "##") |>
    select(subject_id = SAMPLE_NAME,
           population = POPULATION) |>
    distinct()

### Colors
afr_cols <- brewer.pal("Oranges", n = 9)[3:9] |>
    setNames(c("ACB", "ESN", "GWD", "LWK", "MSL", "YRI", "ASW"))

eur_cols <- brewer.pal("Blues", n = 9)[4:8] |>
    setNames(c("GBR", "IBS", "TSI", "CEU", "FIN"))

sas_cols <- brewer.pal("Greens", n = 9)[5:9] |>
    setNames(c("BEB", "PJL", "GIH", "ITU", "STU"))

eas_cols <- brewer.pal("Purples", n = 9)[5:9] |>
    setNames(c("CHB", "CHS", "CDX", "KHV", "JPT"))

amr_cols <- c("lightpink1", "hotpink", "hotpink3", "deeppink") |>
    setNames(c("MXL", "CLM", "PEL", "PUR"))

mgb_cols <- c("black", "grey") |>
    setNames(c("SLE", "Control")) 

all_cols <- c(mgb_cols, afr_cols, eur_cols, sas_cols, eas_cols, amr_cols)

pca_genos <- 
    "/temp_work/ch229163/vcf/prs/allchr.merged.pruned.plink.pca.eigenvec" |>
    read_table() |>
    select(sample_id = 1, PC1:PC10)

pca_mgb <- pca_genos |>
    filter(!grepl("^HG|^NA", sample_id)) |>
    extract(sample_id, c("subject_id"), "[^-]+-(\\d+)", remove = FALSE) |>
    select(subject_id, PC1:PC10)
    
pca_kgp <- pca_genos |>
    filter(grepl("^HG|^NA", sample_id)) |>
    separate(sample_id, c("sample_id", "subject_id"), sep = "_") |>
    select(subject_id, PC1:PC10)

pca_meta <- sle_data |>
    select(subject_id, population = group) |>
    bind_rows(sample_annotation)

pca_df <- 
    bind_rows(pca_mgb, pca_kgp) |>
    inner_join(pca_meta) |>
    mutate(population = factor(population, levels = names(all_cols)))

## UMAP
# Not variance explained in the data because only 10 PCs were approximated
"/temp_work/ch229163/vcf/prs/allchr.merged.pruned.plink.pca.eigenval" |>
    read_lines() |>
    {function(x) tibble(eigenval = as.numeric(x))}() |>
    select(eigenval) |>
    mutate(var_exp = eigenval/sum(eigenval) * 100)

pca_legend_1 <- 
    ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(data = filter(pca_df, population %in% c("Control", "SLE")),
	       aes(color = population), shape = 3) +
    scale_color_manual(values = all_cols[c("Control", "SLE")]) +
    theme_bw() +
    theme(legend.margin = margin(b = -6.5, unit = "cm")) +
    labs(color = "Population:") +
    guides(color = guide_legend(override.aes = list(size = 2, stroke = 1.5)))

pca_legend_1 <- get_legend(pca_legend_1)

pca_legend_2 <- 
    ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(data = filter(pca_df, ! population %in% c("Control", "SLE")),
	       aes(color = population), shape = 1) +
    scale_color_manual(values = all_cols[! names(all_cols) %in% c("Control", "SLE")]) +
    theme_bw() +
    theme(legend.margin = margin(t = -6.5, unit = "cm")) +
    labs(color = NULL) + 
    guides(color = guide_legend(override.aes = list(size = 2, stroke = 1.5)))

pca_legend_2 <- get_legend(pca_legend_2)


plot_pca <- function(mapping) { 
    ggplot(pca_df, mapping) +
    geom_point(data = filter(pca_df, population == "Control"),
	       aes(color = population), shape = 3, size = .75) +
    geom_point(data = filter(pca_df, ! population %in% c("Control", "SLE")),
	       aes(color = population), shape = 1, size = .75) +
    geom_point(data = filter(pca_df, population == "SLE"),
	       aes(color = population), shape = 3, size = .75) +
    scale_color_manual(values = all_cols) +
    scale_x_continuous(breaks = scales::pretty_breaks(4),
		       labels = function(x) round(x, 2)) +
    scale_y_continuous(breaks = scales::pretty_breaks(4),
		       labels = function(x) round(x, 2)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  legend.position = "none")
}

pca_p1 <- plot_pca(aes(PC1, PC2))
pca_p2 <- plot_pca(aes(PC2, PC3))
pca_p3 <- plot_pca(aes(PC3, PC4))
pca_p4 <- plot_pca(aes(PC4, PC5))
pca_p5 <- plot_pca(aes(PC5, PC6))
pca_p6 <- plot_pca(aes(PC6, PC7))
pca_p7 <- plot_pca(aes(PC7, PC8))
pca_p8 <- plot_pca(aes(PC8, PC9))
pca_p9 <- plot_pca(aes(PC9, PC10))


pca_p <- 
    plot_grid(plot_grid(pca_p1, pca_p2, pca_p3, pca_p4, pca_p5, pca_p6, pca_p7, pca_p8, pca_p9, ncol = 3),
	      plot_grid(pca_legend_1, NULL, pca_legend_2, ncol = 1, align = "hv"),
	      nrow = 1, rel_widths = c(1, .25)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/pca.png", pca_p, width = 10, dpi = 600)


mgb_umap <- select(pca_df, PC1:PC7) |>
    umap()

umap_df <- bind_cols(select(pca_df, subject_id, population), as.data.frame(mgb_umap)) |>
    as_tibble()

plot_umap <- 
    ggplot(data = umap_df, 
	   aes(V1, V2, color = population)) +
    geom_point(data = filter(umap_df, population == "Control"),
	       size = .75, shape = 3) +
    geom_point(data = filter(umap_df, ! population %in% c("Control", "SLE")),
	       size = .75, shape = 1) +
    geom_point(data = filter(umap_df, population == "SLE"),
	       size = .75, shape = 3) +
    scale_color_manual(values = all_cols) +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  legend.position = "none") +
    labs(x = "UMAP 1", y = "UMAP 2")

race_colors <- 
    c("White" = eur_cols[[2]], 
      "Spanish" = "cornflowerblue", 
      "Other" = "black",
      "American Indian or Alaska Native" = "pink", 
      "Asian" = "purple", 
      "Black" = "orange",
      "Native Hawaiian or Other Pacific Islander" = "forestgreen",
      "Hispanic or Latino" = "magenta")


plot_umap_race <- umap_df |>
    inner_join(sle_data) |>
    sample_frac(1L) |>
    mutate(subject_id = fct_inorder(subject_id)) |>
    select(subject_id, race, V1, V2) |>
    ggplot(aes(V1, V2, color = race)) +
    geom_point(size = .75, shape = 1) +
    scale_color_manual(values = race_colors,
		       labels = function(x) str_wrap(x, width = 16),
		       na.value = "grey") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(override.aes = list(size = 2, stroke = 1.5)))

ggsave("./plots/umap.png", plot_umap, width = 5, height = 5)
ggsave("./plots/umap_race.png", plot_umap_race, width = 7, height = 5)



