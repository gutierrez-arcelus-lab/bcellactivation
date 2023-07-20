library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(patchwork)

temp_dir <- system("echo $TEMP_WORK", intern = TRUE)

## PCA on genotype data
### Metadata for 1000G samples
sample_meta <- read_tsv("./data/kgp_metadata.tsv")

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

all_cols <- c(afr_cols, eur_cols, sas_cols, eas_cols, amr_cols)


### PCA results

pca_genos <- 
    file.path(temp_dir, "vcf/prs/allchr.1000G.plink.pca.eigenvec") |>
    read_table() |>
    select(sample_name = 1, PC1:PC10) |>
    extract(sample_name, "sample_name", "([^_]+)")

pca_df <- pca_genos |> 
    inner_join(sample_meta, join_by(sample_name)) |>
    mutate(population_code = factor(population_code, levels = names(all_cols))) |>
    select(sample_name, super_population, population_code, everything())

p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = population_code)) +
    geom_point() +
    scale_color_manual(values = all_cols) +
    theme_bw()

p2 <- ggplot(pca_df, aes(x = PC2, y = PC3, color = population_code)) +
    geom_point() +
    scale_color_manual(values = all_cols) +
    theme_bw()

p1o <- p1 + theme(legend.position = "none")
p2o <- p2 + theme(legend.position = "none")

p <- p1o + p2o + get_legend(p1) + plot_layout(nrow = 1, widths = c(1, 1, .33))

ggsave("./pca.png", p, height = 4, width = 10)

selected_amr <- pca_df |>
    top_n(125, -PC3) |>
    select(1:3)

set.seed(1)
ref_panel <- sample_meta |>
    filter(super_population != "AMR", 
	   !population_code %in% c("ACB", "ASW")) |>
    group_by(super_population) |>
    sample_n(125) |>
    ungroup() |>
    bind_rows(selected_amr) |>
    arrange(super_population, population_code, sample_name)

ref_panel |>
    mutate(ids = paste(sample_name, sample_name, sep = "_")) |>
    arrange(sample_name) |>
    pull(ids) |>
    write_lines("./data/kgp_refpanel_ids.txt")

# plot admixture results
ref_panel_admix <- 
    file.path(temp_dir, "vcf/prs/allchr.refpanel.fam") |>
    read_delim(col_names = FALSE, delim = " ") |>
    select(sample_name = X1)

mgb_admix <- 
    file.path(temp_dir, "vcf/prs/allchr.mgb.fam") |>
    read_delim(col_names = FALSE, delim = " ") |>
    select(sample_name = X2)

ref_q <- 
    "./allchr.refpanel.5.Q" |>
    read_delim(delim = " ", col_names = FALSE) |>
    bind_cols(ref_panel_admix) |>
    left_join(sample_meta, join_by(sample_name)) |>
    select(sample_name, population_code, super_population, X1:X5) |>
    pivot_longer(X1:X5, names_to = "ancestry") |>
    mutate(ancestry = fct_inorder(ancestry))

continental_colors <- 
    c("EUR" = all_cols[["GBR"]],
      "AMR" = all_cols[["MXL"]],
      "AFR" = all_cols[["LWK"]],
      "SAS" = all_cols[["PJL"]],
      "EAS" = all_cols[["CHB"]])

ancestry_pop_map <- 
    ref_q |>
    group_by(ancestry, super_population) |>
    summarise(q = mean(value)) |>
    group_by(ancestry) |>
    slice_max(q) |>
    ungroup() |>
    select(ancestry, super_population)

mgb_q <- 
    "./allchr.mgb.5.Q" |>
    read_delim(delim = " ", col_names = FALSE) |>
    bind_cols(mgb_admix) |>
    select(sample_name, X1:X5) |>
    pivot_longer(X1:X5, names_to = "ancestry") |>
    mutate(ancestry = fct_inorder(ancestry)) |>
    left_join(ancestry_pop_map, join_by(ancestry))

sample_order <- mgb_q |>
    group_by(sample_name) |>
    slice_max(value) |>
    ungroup() |>
    arrange(ancestry, value) |>
    pull(sample_name)

admix_p <- 
    ggplot(mgb_q |> mutate(sample_name = factor(sample_name, levels = sample_order)), 
	   aes(x = sample_name, y = value, fill = super_population)) +
    geom_col(width = 1.001) +
    scale_y_continuous(labels = scales::percent,
		       breaks = c(0, .5, 1)) +
    scale_fill_manual(values = continental_colors) +
    facet_grid(.~ancestry, space = "free") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  strip.text = element_blank(),
	  legend.position = "top",
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = NULL, y = "Ancestry percentage")

ggsave("./admix.png", admix_p, width = 12, height = 4)

