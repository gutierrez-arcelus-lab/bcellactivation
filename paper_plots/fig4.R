library(tidyverse)
library(ggrepel)
library(scico)
library(cowplot)
library(glue)
library(tidytext)

# colors
stim_colors <- 
    read_tsv("../figure_colors.txt", col_names = c("stim", "color")) |>
    mutate(stim = str_replace(stim, "_", " "),
	   stim = paste0(stim, "h")) |>
    deframe()

# meta data
donor_ids <- 
    "../atacseq/samplesheet.csv" |>
    read_csv() |>
    mutate(donor_id = basename(fastq_1)) |>
    extract(donor_id, "donor_id", "[^_]+_([^_]+)_.+") |>
    mutate(replicate = paste0("REP", replicate),
	   stim = str_replace(sample, "_", " "),
	   stim = str_replace(stim, "unst", "Unstim"),
	   stim = paste0(stim, "h")) |>
    select(stim, donor_id, donor = replicate) |>
    distinct() |>
    filter(donor_id != "3donors")

# Figure A ####################################################################
pca_file <- "../atacseq/results_deseq2/pcadata_5000peaks.rds"

pca_5000_df <- 
    read_rds(pca_file) |>
    mutate(stim = str_replace(condition, "_", " "),
	   stim = str_replace(stim, "unst", "Unstim"),
	   stim = paste0(stim, "h")) |>
    extract(name, "donor", ".+(REP\\d)") |>
    select(stim, donor, PC1, PC2) |>
    left_join(donor_ids) |>
    mutate(stim = factor(stim, levels = names(stim_colors)))

perc_var <-
    read_rds(pca_file) |>
    attr("percentVar") |>
    {function(x) round(x * 100)}()

fig_a <- 
    ggplot(pca_5000_df, aes(PC1, PC2)) +
    geom_point(aes(fill = stim), shape = 21, size = 3, stroke = .5) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(
	  axis.text = element_text(size = 8),
	  axis.title = element_text(size = 9),
	  legend.text = element_text(size = 8, margin = margin(l = 0)),
	  legend.title = element_text(size = 9),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  panel.grid.major.x = element_line(color = "grey96"),
	  panel.grid.major.y = element_line(color = "grey96"),
	  legend.key.height = unit(.2, "in"), 
	  legend.box.spacing = unit(1, "pt")
	  ) +
    labs(x = sprintf("PC1: %s%% variance", perc_var[1]),
	 y = sprintf("PC2: %s%% variance", perc_var[2]),
	 fill = "Stim:")

fig_a_title <- 
    ggdraw() + 
    draw_label(
	       "PCA shows separation of conditions",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_a_grid <- plot_grid(fig_a_title, fig_a, ncol = 1, rel_heights = c(.1, 1))




# Figure B ####################################################################
da_files <- 
    list.files("../atacseq/results_deseq2",
	       pattern = ".+vs.+\\.tsv$",
	       full.names = TRUE)

da_files <- setNames(da_files, basename(da_files) |> str_remove("\\.tsv"))

da_data <-
    map_dfr(da_files, read_tsv, .id = "comparison") |>
    filter(!is.na(padj)) |>
    separate(comparison, c("stim1", "stim2"), sep = "vs")

da_positive <- 
    da_data |>
    filter(log2FoldChange > 0, padj <= 0.01)

da_negative <- 
    da_data |>
    filter(log2FoldChange < 0, padj <= 0.01)

da_summ <- 
    bind_rows("positive" = da_positive, "negative" = da_negative, .id = "comparison") |>
    mutate(comparison = fct_inorder(comparison)) |>
    count(comparison, stim1, stim2) |>
    mutate_at(vars(stim1, stim2), 
	      ~str_replace(., "_", " ") |> 
	      str_replace("unst", "Unstim") |>
	      paste0("h") |> 
	      factor(levels = names(stim_colors))) |>
   arrange(comparison, stim2, stim1) |>
   mutate(lab = n,
	  n = case_when(comparison == "positive" ~ n,
			comparison == "negative" ~ n * -1),
	  stim_a = case_when(comparison == "positive" ~ stim1,
			     comparison == "negative" ~ stim2),
	  stim_b = case_when(comparison == "positive" ~ stim2,
			     comparison == "negative" ~ stim1)) |>
    select(comparison, stim_a, stim_b, n, lab)

da_plot <- 
    ggplot(da_summ, aes(x = stim_b, y = stim_a)) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = scales::comma(lab)), 
	      color = "white", fontface = "bold", size = 7, size.unit = "pt") +
    scale_x_discrete(position = "top") +
    scale_fill_scico(palette = "bamO", na.value = "white") +
    theme_minimal() + 
    theme(
	  panel.grid = element_blank(),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.ticks.y = element_blank(),
	  legend.position = "none",
	  plot.margin = margin(0, 0, 0, 0, "cm")
    )

cond_vertical <- 
    ggplot(da_summ |> distinct(stim_a),
	   aes(x = factor(1), y = stim_a)) +
    geom_tile(aes(fill = stim_a), color = "black") +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.ticks.y = element_blank(),
	  plot.margin = margin(0, 0, 0, 0)) +
    guides(fill = "none") +
    coord_cartesian(xlim = c(1, 1))

cond_horizontal <- 
    ggplot(da_summ |> distinct(stim_b),
	   aes(y = factor(1), x = stim_b)) +
    geom_tile(aes(fill = stim_b), color = "black") +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.ticks.y = element_blank(),
	  plot.margin = margin(0, 0, 0, 0, "lines")) +
    guides(fill = "none") +
    coord_cartesian(ylim = c(1, 1))


fig_b <- 
    plot_grid(NULL, cond_horizontal, cond_vertical, da_plot,
	      ncol = 2, nrow = 2, align = "hv", rel_heights = c(.1, 1), rel_widths = c(0.05, 1))


fig_b_title <- 
    ggdraw() + 
    draw_label(
	       "# of differentially accessible (DA) peaks",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_b_grid <- 
    plot_grid(fig_b_title, fig_b, ncol = 1, rel_heights = c(.1, 1)) +
    theme(plot.margin = margin(r = 5.5 * 2, unit = "pt"))
	


# Figure C ####################################################################
stims <- read_lines("../atacseq/homer/stims.txt")

homer_df <- 
    glue("../atacseq/homer/results/{stims}/knownResults.txt") |>
    setNames(stims) |>
    map_dfr(~read_tsv(.) |> 
	    janitor::clean_names() |>
	    {function(x) setNames(x, str_remove(names(x), "_of_\\d+$"))}(),
	    .id = "stim") |>
    mutate_at(vars(starts_with("percent_of_")), parse_number) |>
    mutate(stim = factor(stim, levels = stims),
	   fc = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) |>
    filter(!grepl("bias", motif_name, ignore.case = TRUE)) |>
    mutate(motif_name = str_replace(motif_name, "NR/Ini-like", "NR;Ini-like")) |>
    separate(motif_name, c("motif", "dataset", "db"), sep = "/") |>
    separate(motif, c("motif", "tf_family"), sep = "\\(") |>
    mutate(tf_family = str_remove(tf_family, "\\)$")) |>
    mutate(log10p = log_p_value/log(10)) |>
    select(stim, motif, tf_family, dataset, 
	   pct_target =  percent_of_target_sequences_with_motif,
	   pct_bg = percent_of_background_sequences_with_motif,
	   fc, log10p, q_value = q_value_benjamini)

homer_top25 <- 
    homer_df |>
    filter(q_value <= 0.05) |>
    group_by(stim) |>
    top_n(25, -log10p) |>
    ungroup()


homer_top25_inall <- 
    homer_df |>
    inner_join(distinct(homer_top25, motif, tf_family)) |>
    arrange(log10p) |>
    mutate(motif = fct_inorder(motif),
	   motif = fct_rev(motif))

fig_c <- 
    ggplot(homer_top25_inall, 
	   aes(x = stim, y = motif)) +
    geom_point(aes(fill = log2(fc), size = -log10p), 
	       shape = 21, stroke = .2) +
    geom_text(data = filter(homer_top25_inall, q_value <= 0.01), 
	      aes(x = stim, y = motif, label = "*"), 
	      size = 8, fontface = "bold", size.unit = "pt", 
	      nudge_x = .35, nudge_y = -.2) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "left") +
    scale_fill_scico(palette = "vik", midpoint = 0) +
    scale_size(range = c(2, 5)) +
    facet_grid(rows = vars(tf_family), scales = "free", space = "free", switch = "x") +
    theme_minimal() +
    theme(
	  axis.text.x = element_text(size = 8), 
	  axis.text.y = element_text(size = 7), 
	  axis.title = element_blank(),
	  strip.text = element_blank(),
	  strip.clip = "off",
	  panel.spacing.y = unit(.2, "lines"),
	  legend.title = element_text(size = 8),
	  legend.text = element_text(size = 8),
	  legend.box.spacing = unit(-.2, "lines"),
	  ) +
    guides(fill = guide_colorbar(barwidth = .5, barheight = 6)) +
    labs(fill = "log2(FC)", size = "-log10(P)") +
    coord_cartesian(clip = "off")

fig_c_title <- 
    ggdraw() + 
    draw_label(
	       "Motif enrichment of DA peaks in a\ncondition vs. Unstim 24 hours",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_c_grid <- plot_grid(fig_c_title, fig_c, ncol = 1, rel_heights = c(.05, 1))


#homer_heat <- 
#    ggplot(homer_top25_inall, 
#	   aes(x = motif, y = stim)) +
#    geom_tile(aes(fill = log2(fc))) +
#    geom_text(data = filter(homer_top25_inall, q_value <= 0.01), 
#	      aes(x = motif, y = stim, label = "*"), 
#	      size = 7, size.unit = "pt", nudge_y = -.1) +
#    scale_y_discrete(position = "left") +
#    scale_fill_scico(palette = "vik", midpoint = 0) +
#    facet_grid(cols = vars(tf_family), scales = "free", space = "free", switch = "x") +
#    theme_minimal() +
#    theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1, 
#				     margin = margin(0, 0, 0, 0)),
#	  axis.text.y = element_text(size = 8),
#	  axis.title = element_blank(),
#	  strip.text.x.bottom = element_blank(),
#	  panel.spacing.x = unit(0.1, "pt"),
#	  legend.position = "right",
#	  legend.text = element_text(size = 8),
#	  legend.title = element_text(size = 8),
#	  legend.box.spacing = unit(-.25, "lines"),
#	  plot.margin = margin(0, 0, 0, 0),
#	  plot.background = element_rect(fill = "white", color = "white")) +
#    guides(fill = guide_colorbar(barwidth = .25, barheight = 1.5)) +
#    labs(fill = NULL)
#


#homer_plot <- 
#    ggplot(homer_top25, 
#	   aes(x = -log10p, 
#	       y = reorder_within(motif, by = -log10p, within = stim))) +
#    geom_point(aes(fill = log2(fc)), shape = 21, size = 2, stroke = .2) +
#    scale_x_continuous(limits = function(x) c(0, max(x)),
#		       expand = expansion(mult = .1),
#		       breaks = function(x) c(0, max(x)),
#		       labels = function(x) round(x)) +
#    scale_y_reordered() +
#    scale_fill_gradient(low = "white", high = "black",
#			limits = function(x) c(0, max(x)),
#			breaks = function(x) seq(0, max(x), length.out = 3),
#			labels = function(x) round(x)) +
#    facet_wrap(~stim, nrow = 1, scales = "free") +
#    theme_minimal() +
#    theme(axis.text = element_text(size = 8),
#	  axis.title.x = element_text(size = 9, margin = margin(t = -.5, unit = "lines")),
#	  axis.title.y = element_blank(),
#	  strip.text = element_text(size = 9),
#	  strip.clip = "off",
#	  legend.box.spacing = unit(.1, "lines"),
#	  legend.text = element_text(size = 8),
#	  legend.title = element_text(size = 8),
#	  panel.grid.minor = element_blank(),
#	  plot.margin = margin(0, 5.5, 0, 5.5, "pt")
#	  ) +
#    guides(fill = guide_colorbar(barwidth = .25))



# Figure D ####################################################################
traits <- 
    read_tsv("../atacseq/ldsc/data/traits.txt", col_names = c("directory", "trait", "gwas", "ref")) |>
    select(gwas, trait)

results <- 
    read_tsv("../atacseq/ldsc/compiled_results.tsv") |>
    mutate(set = factor(set, levels = c("IL4", "TLR7", "BCR", "DN2")),
	   group = case_when(grepl("T2D|Height|HEIGHTz|Type_2_Diabetes|Schizophrenia|MDD|LDL|HDL|Covid19_Infection|Vaccination|cancer", gwas) ~ "control",
			     TRUE ~ "test")) |>
    left_join(traits, join_by(gwas)) |>
    arrange(trait, gwas, set) |>
    mutate(trait = fct_inorder(trait), gwas = fct_inorder(gwas))


results <- 
    results |> 
    mutate(gwas = case_when(gwas == "GBMI.Asthma" ~ "Asthma (GBMI)",
			    gwas == "UKB_460K.disease_ASTHMA_DIAGNOSED" ~ "Asthma (UKB)",
			    gwas == "PASS_AdultOnsetAsthma_Ferreira2019" ~ "AOA (Ferreira 2019)",
			    gwas == "PASS_ChildOnsetAsthma_Ferreira2019" ~ "COA (Ferreira 2019)",
			    gwas == "UKB_460K.cancer_ALL" ~ "Cancer-ALL (UKB)",
			    gwas == "PASS.Celiac.Dubois2010" ~ "Celiac (Dubois 2010) v4",
			    gwas == "PASS_Celiac" ~ "Celiac (Dubois 2010) v1",
			    gwas == "PASS.Covid19_Infection.hg_v7" ~ "Covid19 infection (hg_v7)",
			    gwas == "PASS_CD_deLange2017" ~ "CD (deLange 2017)",
			    gwas == "PASS_Crohns_Disease" ~ "CD (Jostins 2012)",
			    gwas == "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED" ~ "Allergy/Eczema (UKB)",
			    gwas == "PASS_HDL" ~ "HDL (Teslovich 2010)",
			    gwas == "PASS_Height1" ~ "Height (Lango Allen 2010)",
			    gwas == "UKB_460K.body_HEIGHTz" ~ "Height (UKB)",
			    gwas == "PASS.Height.Yengo2022" ~ "Height (Yengo 2022)",
			    gwas == "PASS_IBD" ~ "IBD (Jostins 2012)",
			    gwas == "PASS_IBD_deLange2017" ~ "IBD (deLange 2017)",
			    gwas == "PASS_LDL" ~ "LDL (Teslovich 2010)",
			    gwas == "PASS_MDD_Howard2019" ~ "MDD (Howard 2019)",
			    gwas == "PASS_MDD_Wray2018" ~ "MDD (Wray 2018)",
			    gwas == "PASS_Multiple_sclerosis" ~ "MS (IMS)",
			    gwas == "PASS_Primary_biliary_cirrhosis" ~ "PBC (Cordell 2015) v1",
			    gwas == "PASS.Primary_Biliary_Cirrhosis.Cordell2015" ~ "PBC (Cordell 2015) v4",
			    gwas == "UKB_460K.disease_PSORIASIS" ~ "Psoriasis (UKB)",
			    gwas == "PASS_Rheumatoid_Arthritis" ~ "RA (Okada 2014)",
			    gwas == "PASS.Rheumatoid_Arthritis.Ishigaki2022" ~ "RA (Ishigaki 2022)",
			    gwas == "PASS.Rheumatoid_Arthritis.Saevarsdottir2022" ~ "RA (Saevarsdottir 2022)",
			    gwas == "PASS.Covid19_Vaccination.Hartonen2023" ~ "Covid19 Vac. (Hartonen 2023)",
			    gwas == "PASS_Schizophrenia" ~ "SCZ (SCZ Consort.)",
			    gwas == "PASS_Schizophrenia_Pardinas2018" ~ "SCZ (Pardinas 2018)",
			    gwas == "PASS_Schizophrenia_Ruderfer2018" ~ "SCZ (Ruderfer 2018)",
			    gwas == "PASS.Schizophrenia.Trubetskoy2022" ~ "SCZ (Trubetskoy 2022)",
			    gwas == "PASS_Lupus" ~ "SLE (Bentham 2015)",
			    gwas == "PASS_Type_1_Diabetes" ~ "T1D (Bradfield 2011)",
			    gwas == "PASS.Type_1_Diabetes.Chiou2021" ~ "T1D (Chiou 2021)",
			    gwas == "PASS_Type_2_Diabetes" ~ "T2D (Morris 2012)",
			    gwas == "UKB_460K.disease_T2D" ~ "T2D (UKB)",
			    gwas == "PASS.Type_2_Diabetes.Xue2018" ~ "T2D (Xue 2018)",
			    gwas == "PASS_Ulcerative_Colitis" ~ "UC (Jostins 2012)",
			    gwas == "PASS_UC_deLange2017" ~ "UC (deLange 2017)"))

fig_d <- 
    ggplot(data = results, aes(x = set, y = gwas)) +
    geom_tile(aes(fill = tau_star)) +
#    geom_text(data = filter(results, pfdr >= 0.01, pfdr <= 0.05),
#	      aes(x = set, y = gwas, label = "*"), 
#	      size = 10, fontface = "bold", size.unit = "pt", nudge_y = -.2) +
    geom_text(data = filter(results, pfdr <= 0.01),
	      aes(x = set, y = gwas, label = "*"), 
	      size = 10, fontface = "bold", size.unit = "pt", nudge_y = -.2) +
    scale_x_discrete(position = "top") +
    scale_fill_scico(palette = "vik", midpoint = 0) +
    facet_grid(rows = vars(group), scales = "free_y", space = "free") +
    theme_minimal() +
    theme(axis.title = element_blank(),
	  axis.text = element_text(size = 8),
	  strip.text = element_blank(),
	  legend.box.spacing = unit(-.25, "lines"),
	  ) +
    labs(fill = "Tau*") +
    guides(fill = guide_colorbar(barheight = 8, barwidth = .5))

fig_d_title <- 
    ggdraw() + 
    draw_label(
	       "DA peaks in stimulated B cells are enriched with\nheritability of immune-mediated diseases",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_d_grid <- plot_grid(fig_d_title, fig_d, ncol = 1, rel_heights = c(.05, 1))

# Final panel
top_grid <- 
    plot_grid(fig_a_grid, NULL, fig_b_grid, nrow = 1, rel_widths = c(1, .05, .85),
	      labels = c("a", "", "b"), label_size = 12)

bottom_grid <- 
    plot_grid(fig_c_grid, NULL, fig_d_grid, nrow = 1, rel_widths = c(.8, .1, 1),
	      labels = c("c", "", "d"), label_size = 12)


final_grid <- 
    plot_grid(top_grid, NULL, bottom_grid, ncol = 1, rel_heights = c(.3, .05, 1)) +
    theme(panel.background = element_rect(color = "white", fill = "white"))

ggsave("fig4.png", final_grid, width = 6.5, height = 8.5)

