library(tidyverse)
library(tidytext)
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

p3 <- ggplot(pca_df, aes(x = PC3, y = PC4, color = population_code)) +
    geom_point() +
    scale_color_manual(values = all_cols) +
    theme_bw()

p4 <- ggplot(pca_df, aes(x = PC4, y = PC5, color = population_code)) +
    geom_point() +
    scale_color_manual(values = all_cols) +
    theme_bw()

p1o <- p1 + theme(legend.position = "none")
p2o <- p2 + theme(legend.position = "none")
p3o <- p3 + theme(legend.position = "none")
p4o <- p4 + theme(legend.position = "none")

ptemp <- p1o + p2o + p3o + p4o + plot_layout(nrow = 2)

p <- plot_grid(ptemp, get_legend(p1), nrow = 1, rel_widths = c(1, .25)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./pca.png", p, height = 7, width = 9)

selected_amr <- pca_df |>
    top_n(125, -PC3) |>
    select(1:3)

selected_afr <- pca_df |>
    top_n(125, -PC1) |>
    select(1:3)

selected_sas <- pca_df |>
    top_n(125, PC3) |>
    select(1:3)


set.seed(1)
ref_panel <- sample_meta |>
    filter(super_population %in% c("EAS", "EUR")) |>
    group_by(super_population) |>
    sample_n(125) |>
    ungroup() |>
    bind_rows(selected_amr, selected_afr, selected_sas) |>
    arrange(super_population, population_code, sample_name)

ref_panel |>
    mutate(ids = paste(sample_name, sample_name, sep = "_")) |>
    arrange(sample_name) |>
    pull(ids) |>
    write_lines("./data/kgp_refpanel_ids.txt")

## Run ADMIXTURE


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
    pivot_longer(X1:X5, names_to = "k") |>
    mutate(k = fct_inorder(k))

continental_colors <- 
    c("EUR" = all_cols[["GBR"]],
      "AMR" = all_cols[["MXL"]],
      "AFR" = all_cols[["MSL"]],
      "SAS" = all_cols[["PJL"]],
      "EAS" = all_cols[["CDX"]])

ancestry_pop_map <- 
    ref_q |>
    group_by(k, super_population) |>
    summarise(q = mean(value)) |>
    group_by(k) |>
    slice_max(q) |>
    ungroup() |>
    select(k, ancestry = super_population)

mgb_q <- 
    "./allchr.mgb.5.Q" |>
    read_delim(delim = " ", col_names = FALSE) |>
    bind_cols(mgb_admix) |>
    select(sample_name, X1:X5) |>
    pivot_longer(X1:X5, names_to = "k") |>
    mutate(k = fct_inorder(k)) |>
    left_join(ancestry_pop_map, join_by(k)) |>
    select(-k)

top_ancestry <- mgb_q |>
    group_by(sample_name) |>
    slice_max(value) |>
    ungroup() |>
    count(top = ancestry, sort = TRUE)

top_ancestry_within <- 
    mgb_q |>
    group_by(sample_name) |>
    mutate(max_ancestry = ancestry[which.max(value)]) |>
    ungroup() |>
    group_by(max_ancestry, ancestry) |>
    summarise(m = mean(value)) |>
    ungroup() |>
    mutate(max_ancestry = factor(max_ancestry, levels = top_ancestry$top)) |>
    arrange(max_ancestry, desc(m)) |>
    {function(x) split(x, x$max_ancestry)}() |>
    map("ancestry")

tmp <- 
    mgb_q |>
    group_by(sample_name) |>
    mutate(max_ancestry = ancestry[which.max(value)]) |>
    ungroup() |>
    mutate(max_ancestry = factor(max_ancestry, levels = top_ancestry$top)) |>
    {function(x) split(x, x$max_ancestry)}()

plot_df <- 
    map_dfr(top_ancestry$top, 
	   function(x) {
	    
	    ancestries <- top_ancestry_within[[x]] 

	    tmp[[x]] |>
	    pivot_wider(names_from = ancestry, values_from = value) |>
	    select(sample_name, max_ancestry, !!!ancestries) |>
	    arrange(desc(.data[[ancestries[1]]]),
		    .data[[ancestries[2]]],
		    .data[[ancestries[3]]],
		    .data[[ancestries[4]]],
		    .data[[ancestries[5]]])
	   }) |>
    mutate(sample_name = fct_inorder(sample_name)) |>
    pivot_longer(-c(1:2), names_to = "ancestry") |>
    mutate_at(vars(max_ancestry, ancestry), ~factor(., levels = top_ancestry$top)) |>
    arrange(sample_name, ancestry)

admix_p <- 
    ggplot(plot_df,
	   aes(x = sample_name, 
	       y = value, 
	       fill = ancestry)) +
    geom_col(width = .80, position = "fill") +
    scale_y_continuous(labels = scales::percent,
		       breaks = c(0, .5, 1),
		       limits = c(0, 1)) +
    scale_fill_manual(values = continental_colors) +
    facet_grid(. ~ max_ancestry, space = "free_x", scales = "free_x") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  strip.text = element_blank(),
	  legend.position = "top",
	  plot.background = element_rect(fill = "white", color = "white"),
	  panel.background = element_rect(fill = "white", color = "white"),
	  panel.grid = element_blank()) +
    labs(x = NULL, y = "Ancestry percentage", fill = "Ancestry")

ggsave("./admix.png", admix_p, width = 12, height = 4)

# save admixture results
plot_df |>
    select(subject_id = sample_name, ancestry, value) |>
    pivot_wider(names_from = "ancestry", values_from = "value") |>
    mutate(subject_id = sub("([^-]+-)", "", subject_id)) |>
    write_tsv("./admix_results.tsv")

# plot admixture profiles per self-reported race
sle_data <- read_tsv("../mgb_data/sle_data.tsv", col_types = c(subject_id = "c"))

admix_race_df <- plot_df |>
    mutate(sample_name = as.character(sample_name),
	   subject_id = sub("^[^-]+-(\\d+)$", "\\1", sample_name)) |>
    left_join(sle_data, join_by(subject_id)) |>
    select(subject_id, race, group, ancestry, q = value) |>
    drop_na()

race_order <- admix_race_df |> 
    count(race, sort = TRUE) |>
    pull(race)

admix_race_df <- admix_race_df |>
    mutate(race = factor(race, levels = race_order)) |>
    arrange(race, ancestry, )

admix_race_df |> filter(race == "Hispanic or Latino") |> count(group)
admix_race_df |> filter(group == "SLE", ancestry == "AMR") |> arrange(desc(q))
