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

# Not variance explained in the data because only 10 PCs were approximated
#"/temp_work/ch229163/vcf/prs/allchr.merged.pruned.plink.pca.eigenval" |>
#    read_lines() |>
#    {function(x) tibble(eigenval = as.numeric(x))}() |>
#    select(eigenval) |>
#    mutate(var_exp = eigenval/sum(eigenval) * 100)

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
    theme(legend.margin = margin(t = -6.5, unit = "cm"),
	  legend.spacing.y = unit(0, 'cm')) +
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

# UMAP
continental_colors <- 
    c("AFR" = all_cols[["MSL"]],
      "AMR" = all_cols[["PEL"]],
      "EAS" = all_cols[["CHS"]],
      "EUR" = all_cols[["GBR"]],
      "SAS" = all_cols[["BEB"]])

race_colors <- 
    c("White" = "cornflowerblue", 
      "Spanish" = "midnightblue", 
      "Other" = "grey",
      "American Indian or Alaska Native" = "black", 
      "Asian" = "purple", 
      "Black" = "brown",
      "Native Hawaiian or Other Pacific Islander" = "goldenrod",
      "Hispanic or Latino" = "hotpink")

shapes <- 
    tibble(group = levels(umap_df$population)) |>
    mutate(shape = case_when(group %in% c("Control", "SLE") ~ 3L,
			     TRUE ~ 1L)) |>
    deframe()

set.seed(1)
mgb_umap <- select(pca_df, PC1:PC8) |>
    umap(n_threads = future::availableCores())

umap_df <- 
    bind_cols(select(pca_df, subject_id, population), 
	      as.data.frame(mgb_umap)) |>
    as_tibble()
    
umap_df <- umap_df |>
    mutate(population = fct_relevel(population, "SLE", after = Inf)) |>
    arrange(population) |>
    mutate(subject_id = fct_inorder(subject_id))

legend_order <- 
    umap_df$population |>
    unique() |>
    fct_inorder() |>
    fct_relevel("SLE", after = 1L) |>
    sort() |>
    as.character()

plot_umap <- 
    ggplot(data = umap_df,
	   aes(V1, V2, color = population, shape = population)) +
	geom_point(size = 1) +
	scale_color_manual(values = all_cols, 
			   breaks = legend_order) +
	scale_shape_manual(values = shapes,
			   breaks = legend_order) +
	theme_minimal() +
	theme(panel.grid = element_blank(),
	      axis.text = element_blank(),
	      plot.title = element_text(size = 12),
	      plot.margin = margin(.5, .5, 1, .5, unit = 'cm'),
	      plot.background = element_rect(color = "white", fill = "white")) +
	labs(x = "UMAP 1", 
	     y = "UMAP 2", 
	     color = "Group:", 
	     shape = "Group:",
	     title = "UMAP with 8 PCs using\n1000 Genomes as reference panel") +
	guides(color = guide_legend(override.aes = list(size = 2, stroke = 1.25)))


umap_race_df <- umap_df |>
    inner_join(sle_data) |>
    sample_frac(1L) |>
    mutate(subject_id = fct_inorder(subject_id)) |>
    select(subject_id, race, V1, V2) 

plot_umap_race <-
    ggplot(umap_race_df, 
	   aes(V1, V2, color = race)) +
    geom_point(size = 1, shape = 1) +
    scale_color_manual(values = race_colors,
		       labels = function(x) str_wrap(x, width = 16),
		       na.value = "grey") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  axis.text = element_blank(),
	  plot.title = element_text(size = 12),
	  plot.margin = margin(.5, .5, 1, .5, unit = 'cm'),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = "UMAP 1", 
	 y = "UMAP 2", 
	 title = "UMAP colored by reported race\nin MGB Biobank") +
    guides(color = guide_legend(override.aes = list(size = 2, stroke = 1.25)))

# Ancestry proportions from ADMIXTURE
admix <- 
    "./admixture/admix_results.tsv" |>
    read_tsv(col_types = c(subject_id = "c")) |>
    left_join(umap_prs_df) |>
    select(subject_id, group, race, EUR:SAS, prs) |>
    pivot_longer(EUR:SAS, names_to = "ancestry")
    

## Check UMAP colored by ancestry proportions
umap_anc_df <- 
    filter(umap_df, population %in% c("SLE", "Control")) |>
    left_join(select(admix, subject_id, ancestry, value), join_by(subject_id))

umaps_anc <- 
    map(names(continental_colors), 
	function(anc) {
	    ggplot(data = filter(umap_anc_df, ancestry == anc),
		   aes(x = V1, y = V2)) +
	geom_point(aes(color = value), size = .25) +
	scale_color_gradient(low = "grey96", high = continental_colors[anc],
			     guide = guide_colorbar(barheight = .25),
			     breaks = c(0, .5, 1),
			     limits = c(0, 1)) +
	facet_wrap(~ancestry) +
	theme_bw() +
	theme(axis.text = element_blank(),
	      axis.ticks = element_blank(),
	      panel.grid = element_blank(),
	      legend.position = "bottom",
	      legend.margin = margin(t = -.25, b = 0, unit = "cm"),
	      plot.margin = margin(t = .5, unit = "cm"),
	      plot.background = element_rect(color = "white", fill = "white")) +
	labs(x = NULL, y = NULL, color = NULL)}) |>
    {function(x) plot_grid(plotlist = x, nrow = 1)}() +
    labs(title = "Proportion of each ancestry as estimated with ADMIXTURE") +
    theme(plot.title = element_text(size = 12, hjust = 0),
	  plot.background = element_rect(color = "white", fill = "white"))
       
umap_out <- 
    plot_grid(
	      plot_grid(plot_umap, plot_umap_race, nrow = 1, labels = c("A)", "B)")),
	      plot_grid(NULL, umaps_anc, NULL, 
			nrow = 1, rel_widths = c(.1, 1, .1),
			labels = c("C)")),
	      ncol = 1, rel_heights = c(1, .5)) +
    theme(plot.background = element_rect(color = "white", fill = "white"))

ggsave("plots/umaps.png",
       umap_out,
       width = 11, height = 7, dpi = 600)



# PRS
prs_df <- 
    "./common_data/prs_p+t_5e-4.tsv" |>
    read_tsv(col_type = c("subject_id" = "c"))

umap_prs_df <- umap_df |>
    filter(population %in% c("SLE", "Control")) |>
    left_join(prs_df, join_by(subject_id)) 

admix_filt <- admix |>
    filter(value > 0.1)

p1 <- 
    ggplot(data = admix_filt,
	   aes(x = value, y = prs)) +
	geom_point(data = filter(admix_filt, group == "Control"),
		   aes(color = ancestry), size = .5, alpha = .25) +
	geom_point(data = filter(admix_filt, group == "SLE"),
		   aes(fill = ancestry), size = 1.5, alpha = 1, 
		   shape = 21, stroke = .25) +
	geom_smooth(se = FALSE, method = "lm", 
		    color = "black", alpha = .5) +
	scale_x_continuous(limits = c(.1, 1L),
			   breaks = c(.1, .5, 1)) +
	scale_color_manual(values = continental_colors) +
	scale_fill_manual(values = continental_colors) +
	facet_wrap(~ancestry, nrow = 1) +
	theme_bw() +
	theme(panel.grid = element_blank(),
	      legend.position = "none") +
	labs(x = "Proportion of ancestry", y = "PRS")

ggsave("./plots/prs_ancestry.png", p1, height = 2, width = 8)

prs_violin <- admix |>
    filter(value >= .75) |>
    ggplot(aes(x = group, y = prs, color = ancestry)) +
	geom_quasirandom(method = "smiley", alpha = .5, size = .75) +
	scale_color_manual(values = continental_colors) +
	facet_wrap(~ancestry, nrow = 1) +
	theme_bw() +
	theme(legend.position = "none",
	      panel.grid.major.x = element_blank(),
	      panel.grid.minor.y = element_blank(),
	      panel.grid.major.y = element_line(color = "grey96")) +
	labs(x = NULL)
     
ggsave("./plots/prs_violin.png", prs_violin, height = 2, width = 8)

prs_race_df <- admix |>
    drop_na(race)

prs_race <-
    ggplot(data = prs_race_df, aes(x = str_wrap(race, 20), y = prs)) +
	geom_quasirandom(method = "smiley", alpha = .5, size = .75) +
	theme_bw() +
	theme(panel.grid.major.x = element_blank(),
	      panel.grid.minor.y = element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
	labs(x = NULL, y = "PRS")

ggsave("./plots/prs_race.png", prs_race, height = 5, width = 5)

prs_race_df |>
    group_by(race, ancestry) |>
    summarise(m = mean(value)) |>
    ungroup() |>
    arrange(race, desc(m)) |>
    print(n = Inf)


# AUC
library(pROC)

roc_full <- roc(prs_df$group, prs_df$prs)

roc_df <- 
    roc_full[c("sensitivities", "specificities")] |>
    as.data.frame() |> 
    as_tibble() |>
    mutate(ancestry = "Total") |>
    select(ancestry, everything())

rocs_ancestries <- admix |>
    filter(value >= .75) |> 
    filter(ancestry %in% c("AFR", "EAS", "EUR")) |>
    {function(x) split(x, x$ancestry)}() |>
    map(function(x) roc(x$group, x$prs)) |>
    map(smooth) |>
    map(function(x) x[c("sensitivities", "specificities")]) |>
    map(as.data.frame) |>
    map_df(as_tibble, .id = "ancestry")

roc_colors <- c("Total" = "black", continental_colors[c("AFR", "EAS", "EUR")])

roc_p <- 
    bind_rows(roc_df, rocs_ancestries) |>
    mutate(ancestry = fct_inorder(ancestry)) |>
    ggplot(aes(x = specificities, y = sensitivities, color = ancestry)) +
    geom_line(linetype = 2, linewidth = 1.25) +
    scale_x_reverse() +
    scale_color_manual(values = roc_colors) +
    theme_minimal() +
    theme(panel.grid = element_line(color = "grey96"),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = "Specificity", y = "Sensitivity", color = "Group:")

auc_df <- admix |>
    filter(value >= .75) |> 
    filter(ancestry %in% c("AFR", "EAS", "EUR")) |>
    {function(x) split(x, x$ancestry)}() |>
    map(function(x) roc(x$group, x$prs)) |>
    map(auc) |>
    map(as.numeric) |>
    c("Total" = auc(roc_full)) |>
    bind_rows(.id = "ancestry") |>
    pivot_longer(everything(), names_to = "ancestry", values_to = "auc")

auc_p <- 
    ggplot(data = auc_df, aes(ancestry, auc)) +
    geom_col(aes(fill = ancestry)) +
    scale_fill_manual(values = roc_colors) +
    scale_y_continuous(limits = c(0, 1),
		       breaks = seq(0, 1, .2)) +
    theme_minimal() +
    theme(legend.position = "none",
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(y = "AUC")

ggsave("./plots/roc.png", 
       (roc_p | auc_p) + plot_layout(nrow = 1, widths = c(1, .5)), 
       width = 8, height = 4)



