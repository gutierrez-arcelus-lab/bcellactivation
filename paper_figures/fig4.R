# ==============================================================================
# Description: Generates Figure 4 (ATAC-seq).
#              Panel A: PCA of ATAC-seq samples across conditions.
#              Panel B: Number of differentially accessible (DA) peaks.
#              Panel C: Homer motif enrichment for DA peaks (filtered by RNA expression).
#              Panel D: S-LDSC heritability enrichment of immune/control GWAS traits.
# ==============================================================================

library(tidyverse)
library(ggrepel)
library(scico)
library(cowplot)
library(glue)
library(tidytext)

# Universal helper function to guarantee perfect title/label alignment
# Anchors text to the absolute top (y=1) and pushes it 15pt right to clear the letter label
create_title <- function(text) {
    ggdraw() + 
    draw_label(text, x = 0, y = 1, vjust = 1, hjust = 0, size = 7) + 
    theme(plot.margin = margin(t = 5, b = 7, l = 15, unit = "pt"))
}

# -----------------------------------------------------------------------------
# GLOBAL AESTHETICS & SETTINGS
# -----------------------------------------------------------------------------

# Set global ggplot theme to enforce max 7pt font size
theme_set(theme_minimal(base_size = 7))

# Load color palette and format condition names to match upstream scripts
stim_colors <- 
    read_tsv("./figure_colors_final.txt") |>
    filter(Time %in% c(0, 24),
	   Condition %in% c("Unstim", "IL-4c", "TLR7c", "BCRc", "DN2c")) |>
    unite("stim", c(Condition, Time), sep = " ") |>
    mutate(stim = paste0(stim, "h")) |>
    select(stim, Hex) |>
    deframe()

# Extract metadata from samplesheet using regex
donor_ids <- 
    "../03_atacseq/1-processing/data/samplesheet.csv" |>
    read_csv() |>
    mutate(donor_id = basename(fastq_1)) |>
    extract(donor_id, "donor_id", "[^_]+_([^_]+)_.+") |>
    mutate(replicate = paste0("REP", replicate),
	   stim = str_replace(sample, "_", " "),
	   stim = str_replace(stim, "unst", "Unstim"),
	   stim = str_replace(stim, "IL4", "IL-4c"),
	   stim = str_replace(stim, "BCR", "BCRc"),
	   stim = str_replace(stim, "TLR7", "TLR7c"),
	   stim = str_replace(stim, "DN2", "DN2c"),
	   stim = paste0(stim, "h")) |>
    select(stim, donor_id, donor = replicate) |>
    distinct() |>
    filter(donor_id != "3donors")

# -----------------------------------------------------------------------------
# Fig A: PCA of samples
# -----------------------------------------------------------------------------
pca_file <- 
    "../03_atacseq/2-differential_peaks/results/pcadata_5000peaks.rds"

pca_5000_df <- 
    read_rds(pca_file) |>
    mutate(stim = str_replace(condition, "_", " "),
	   stim = str_replace(stim, "unst", "Unstim"),
	   stim = str_replace(stim, "IL4", "IL-4c"),
	   stim = str_replace(stim, "BCR", "BCRc"),
	   stim = str_replace(stim, "TLR7", "TLR7c"),
	   stim = str_replace(stim, "DN2", "DN2c"),
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
    theme(
	  axis.text = element_text(size = 7),
	  axis.title = element_text(size = 7),
	  legend.text = element_text(size = 7, margin = margin(l = 0)),
	  legend.title = element_text(size = 7),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  panel.grid.major.x = element_line(color = "grey96"),
	  panel.grid.major.y = element_line(color = "grey96"),
	  legend.key.height = unit(.1, "in"), 
	  legend.box.spacing = unit(1, "lines"),
	  plot.margin = margin(0.5, 0, 0, 0.5, "lines")
	  ) +
    coord_cartesian(clip = "off") +
    labs(x = sprintf("PC1 (%s%%)", perc_var[1]),
	 y = sprintf("PC2 (%s%%)", perc_var[2]),
	 fill = "Stim:")

fig_a_title <- create_title("PCA shows separation of conditions")

fig_a_grid <- 
    plot_grid(fig_a_title, fig_a,
	      ncol = 1, rel_heights = c(.1, 1))

# -----------------------------------------------------------------------------
# Fig B: Number of differentially accessible peaks
# -----------------------------------------------------------------------------

# Dynamically load all differential accessibility comparison results
da_files <- 
    list.files("../03_atacseq/2-differential_peaks/results/",
	       pattern = ".+vs.+\\.tsv$",
	       full.names = TRUE)

da_files <- setNames(da_files, basename(da_files) |> str_remove("\\.tsv"))

da_data <-
    map_dfr(da_files, read_tsv, .id = "comparison") |>
    filter(!is.na(padj)) |>
    separate(comparison, c("stim1", "stim2"), sep = "vs")

# Separate into closing and opening peaks based on log2FoldChange
da_positive <- 
    da_data |>
    filter(log2FoldChange > 0, padj <= 0.01)

da_negative <- 
    da_data |>
    filter(log2FoldChange < 0, padj <= 0.01)

# Compile DA counts. Negative peaks are multiplied by -1 to plot on a diverging scale.
da_summ <- 
    bind_rows("positive" = da_positive, "negative" = da_negative, .id = "comparison") |>
    mutate(comparison = fct_inorder(comparison)) |>
    count(comparison, stim1, stim2) |>
    mutate_at(vars(stim1, stim2), 
	      ~str_replace(., "_", " ") |> 
	      str_replace("unst", "Unstim") |>
	      str_replace("IL4", "IL-4c") |>
	      str_replace("BCR", "BCRc") |>
	      str_replace("TLR7", "TLR7c") |>
	      str_replace("DN2", "DN2c") |>
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
	      color = "black", fontface = "bold", size = 7, size.unit = "pt") +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient2(high = '#e66101', mid = 'white', low = '#5e3c99', midpoint = 0) +
    theme(
	  panel.grid = element_blank(),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.ticks.y = element_blank(),
	  legend.position = "none",
	  plot.margin = margin(0, 0, 0, 0, "lines")
    )

cond_vertical <- 
    ggplot(da_summ |> distinct(stim_a),
	   aes(x = factor(1), y = stim_a)) +
    geom_tile(aes(fill = stim_a), color = "black", linewidth = .5) +
    scale_fill_manual(values = stim_colors) +
    theme(panel.grid = element_blank(),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.ticks.y = element_blank(),
	  plot.margin = margin(0, 0.05, 0, 0, "lines")) +
    guides(fill = "none") +
    coord_cartesian(xlim = c(0.5, 1.5), expand = FALSE)

cond_horizontal <- 
    ggplot(da_summ |> distinct(stim_b),
	   aes(y = factor(1), x = stim_b)) +
    geom_tile(aes(fill = stim_b), color = "black", linewidth = .5) +
    scale_fill_manual(values = stim_colors) +
    theme(panel.grid = element_blank(),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.ticks.y = element_blank(),
	  plot.margin = margin(0, 0, 0.05, 0, "lines")) +
    guides(fill = "none") +
    coord_cartesian(ylim = c(0.5, 1.5), expand = FALSE)

fig_b <- 
    plot_grid(NULL, cond_horizontal, cond_vertical, da_plot,
	      ncol = 2, nrow = 2, align = "hv",
	      rel_heights = c(.125, 1), rel_widths = c(0.05, 1)) +
    theme(plot.margin = margin(0.5, 0.5, 0, 0, "lines"))

fig_b_title <- create_title("Number of differentially accessible (DA) peaks")

fig_b_grid <- 
    plot_grid(fig_b_title, fig_b, ncol = 1, rel_heights = c(.1, 1))
	
# -----------------------------------------------------------------------------
# Fig C: Homer motif enrichment for DA peaks
# -----------------------------------------------------------------------------

# Load transcription factor motifs mapped to genes 
motif_genes <- 
    "../03_atacseq/3-motif_analysis/results/homer_motif_genes.tsv" |>
    read_tsv()

# Filter ATAC-seq motifs to ensure the corresponding TF is physically expressed 
# in the matched RNA-seq dataset (Mean TPM >= 5 in at least one condition)
gene_expression <-
    "../01_rnaseq_lowinput/1_quantification/results/expression_pooled_reps.tsv" |>
    read_tsv()

rnaseq_meta <- 
    gene_expression |>
    distinct(sample_id) |>
    separate(sample_id, c("donor", "stim", "timep"), sep = "_", remove = FALSE) |>
    filter(stim %in% c("IL4", "TLR7", "BCR", "DN2"),
	   timep %in% c("4hrs", "24hrs"))

genes_expressed <- 
    gene_expression |>
    inner_join(rnaseq_meta, join_by(sample_id)) |>
    summarize(mean_tpm = mean(tpm), 
	      .by = c(gene_id, gene_name, stim, timep)) |>
    filter(any(mean_tpm >= 5), 
	   .by = c(gene_id, gene_name)) |>
    distinct(gene_id, gene_name) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

selected_motifs <-
    motif_genes |>
    mutate(struct = case_when(grepl(":", gene_symbol) ~ "complex",
			      TRUE ~ "single")) |>
    separate_rows(gene_symbol, sep = ", ") |>
    group_by(motif_name) |>
    mutate(group_id = 1:n()) |>
    ungroup() |>
    separate_rows(gene_symbol, sep = ":") |>
    mutate(is_expressed = gene_symbol %in% genes_expressed$gene_name) |>
    group_by(motif_name, group_id) |>
    filter(all(is_expressed)) |>
    ungroup() |>
    distinct(motif_name)

homer_df <-
    "../03_atacseq/3-motif_analysis/results/results.tsv" |>
    read_tsv()

homer_top <- 
    homer_df |>
    inner_join(selected_motifs) |>
    group_by(tf_name) |>
    slice_min(log10p) |>
    ungroup() |> 
    top_n(50, log10p*-1)

motif_family_order <- 
    homer_top |>
    group_by(tf_family) |>
    slice_min(log10p) |>
    arrange(log10p) |>
    select(motif_name, tf_family, fc, log10p, q_value) |>
    pull(tf_family)

homer_top_inall <- 
    homer_df |>
    inner_join(distinct(homer_top, motif_name, tf_family)) |>
    arrange(log10p) |>
    mutate(stim = 
	   factor(stim, levels = c("IL-4c", "TLR7c", "BCRc", "DN2c")),
	   tf_name = fct_inorder(tf_name),
	   tf_name = fct_rev(tf_name),
	   tf_family = factor(tf_family, levels = motif_family_order))

fig_c <- 
    ggplot(homer_top_inall, 
	   aes(x = stim, y = tf_name)) +
    geom_point(aes(fill = log2(fc), size = -log10p), 
	       shape = 21, stroke = .2) +
    geom_text(data = filter(homer_top_inall, q_value <= 0.01), 
	      aes(x = stim, y = tf_name, label = "*"), 
	      size = 10, fontface = "bold", size.unit = "pt", 
	      nudge_x = .25, vjust = 0.8) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "left") +
    scale_fill_gradient2(high = '#e66101', mid = 'white', low = '#5e3c99', midpoint = 0) +
    scale_size(range = c(2, 5)) +
    facet_grid(rows = vars(tf_family), scales = "free", space = "free", switch = "x") +
    theme(
	  axis.text.x.top = element_text(size = 7, hjust = 0.5), 
	  axis.text.y = element_text(size = 7), 
	  axis.title = element_blank(),
	  strip.text = element_blank(),
	  strip.clip = "off",
	  panel.spacing.y = unit(.2, "lines"),
	  legend.title = element_text(size = 7),
	  legend.text = element_text(size = 7),
	  legend.box.spacing = unit(-.2, "lines"),
	  plot.margin = margin(0.5, 0, 0, 0.5, "lines")
	  ) +
    guides(fill = guide_colorbar(barwidth = .5, barheight = 6)) +
    labs(fill = "log2(FC)", size = "-log10(P)") +
    coord_cartesian(clip = "off")

fig_c_title <- create_title("Motif enrichment of DA peaks in a condition\nvs. Unstim 24 hours")

fig_c_grid <- plot_grid(fig_c_title, fig_c, ncol = 1, rel_heights = c(.05, 1))

# -----------------------------------------------------------------------------
# Fig D: S-LDSC heritability enrichment
# -----------------------------------------------------------------------------
traits <- 
    read_tsv("../03_atacseq/4-heritability/data/traits.txt", col_names = c("directory", "trait", "gwas", "ref")) |>
    select(gwas, trait)

# Group GWAS traits into "control" sets (height, etc.) 
# and "test" sets (immune-mediated diseases)
results <- 
    read_tsv("../03_atacseq/4-heritability/compiled_results.tsv") |>
    rename(gwas = trait) |>
    mutate(set = paste0(set, "c"),
	   set = factor(set, levels = c("IL4c", "TLR7c", "BCRc", "DN2c")),
	   group = case_when(grepl("T2D|Height|HEIGHTz|Type_2_Diabetes|Schizophrenia|MDD|LDL|HDL|cancer", gwas) ~ "control",
			     TRUE ~ "test")) |>
    left_join(traits, join_by(gwas)) |>
    filter(!grepl("cancer_ALL", gwas)) |>
    mutate(trait = fct_inorder(trait), 
	   gwas = fct_inorder(gwas),
	   group = fct_inorder(group))

results <- 
    results |> 
    mutate(gwas = recode(gwas,
			 "GBMI.Asthma" = "Asthma (GBMI)",
			 "UKB_460K.disease_ASTHMA_DIAGNOSED" = "Asthma (UKB)",
			 "PASS_AdultOnsetAsthma_Ferreira2019" = "Asthma-adult (Ferreira2019)",
			 "PASS_ChildOnsetAsthma_Ferreira2019" = "Asthma-child (Ferreira2019)",
			 "PASS_Celiac" = "Celiac (Dubois2010)",
			 "PASS.Covid19_Infection.hg_v7" = "Covid19 infection (hg_v7)",
			 "PASS_CD_deLange2017" = "CD (deLange2017)",
			 "PASS_Crohns_Disease" = "CD (Jostins2012)",
			 "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED" = "Allergy/Eczema (UKB)",
			 "PASS_HDL" = "HDL (Teslovich 2010)",
			 "PASS_Height1" = "Height (Lango Allen2010)",
			 "UKB_460K.body_HEIGHTz" = "Height (UKB)",
			 "PASS.Height.Yengo2022" = "Height (Yengo2022)",
			 "PASS_IBD" = "IBD (Jostins 2012)",
			 "PASS_IBD_deLange2017" = "IBD (deLange2017)",
			 "PASS_LDL" = "LDL (Teslovich 2010)",
			 "PASS_MDD_Howard2019" = "MDD (Howard2019)",
			 "PASS_MDD_Wray2018" = "MDD (Wray 2018)",
			 "PASS_Multiple_sclerosis" = "MS (IMS)",
			 "PASS_Primary_biliary_cirrhosis" = "PBC (Cordell2015)",
			 "UKB_460K.disease_PSORIASIS" = "Psoriasis (UKB)",
			 "PASS_Rheumatoid_Arthritis" = "RA (Okada 2014)",
			 "PASS.Rheumatoid_Arthritis.Ishigaki2022" = "RA (Ishigaki2022)",
			 "PASS.Rheumatoid_Arthritis.Saevarsdottir2022" = "RA (Saevarsdottir2022)",
			 "PASS_Schizophrenia" = "SCZ (SCZ Consort.)",
			 "PASS_Schizophrenia_Pardinas2018" = "SCZ (Pardinas2018)",
			 "PASS_Schizophrenia_Ruderfer2018" = "SCZ (Ruderfer2018)",
			 "PASS.Schizophrenia.Trubetskoy2022" = "SCZ (Trubetskoy2022)",
			 "PASS_Lupus" = "SLE (Bentham2015)",
			 "PASS_Type_1_Diabetes" = "T1D (Bradfield2011)",
			 "PASS.Type_1_Diabetes.Chiou2021" = "T1D (Chiou2021)",
			 "PASS_Type_2_Diabetes" = "T2D (Morris2012)",
			 "UKB_460K.disease_T2D" = "T2D (UKB)",
			 "PASS.Type_2_Diabetes.Xue2018" = "T2D (Xue2018)",
			 "PASS_Ulcerative_Colitis" = "UC (Jostins2012)",
			 "PASS_UC_deLange2017" = "UC (deLange2017)")) |>
    arrange(group, gwas) |>
    mutate(gwas = factor(gwas),
           gwas = fct_relevel(gwas, "Covid19 infection (hg_v7)", after = Inf),
           gwas = fct_rev(gwas))

fig_d <- 
    ggplot(data = results, aes(x = set, y = gwas)) +
    geom_tile(aes(fill = tau_star)) +
    geom_text(data = filter(results, pfdr <= 0.01),
	      aes(x = set, y = gwas, label = "*"), 
	      size = 10, fontface = "bold", size.unit = "pt", nudge_y = -.2) +
    scale_x_discrete(position = "top") +
    scale_fill_gradient2(high = '#e66101', mid = 'white', low = '#5e3c99', midpoint = 0) +
    facet_grid(rows = vars(group), scales = "free_y", space = "free") +
    theme(axis.title = element_blank(),
	  axis.text = element_text(size = 7),
	  strip.text = element_blank(),
	  legend.box.spacing = unit(-.25, "lines"),
	  legend.text = element_text(size = 7),
	  legend.title = element_text(size = 7),
	  plot.margin = margin(0.5, 0.5, 0, 0, "lines")
	  ) +
    labs(fill = "Tau*") +
    guides(fill = guide_colorbar(barheight = 5, barwidth = .5))

fig_d_title <- 
    create_title("DA peaks in stimulated B cells are enriched with\nheritability of immune-mediated diseases")

fig_d_grid <- plot_grid(fig_d_title, fig_d, ncol = 1, rel_heights = c(.05, 1))

# -----------------------------------------------------------------------------
# Compilation & Export
# -----------------------------------------------------------------------------
top_grid <- 
    plot_grid(fig_a_grid, NULL, fig_b_grid, nrow = 1, rel_widths = c(1, .1, 1),
	      labels = c("a", " ", "b"), label_size = 10)

bottom_grid <- 
    plot_grid(fig_c_grid, NULL, fig_d_grid, nrow = 1, rel_widths = c(1, .1, 1),
	      labels = c("c", " ", "d"), label_size = 10)

final_grid <- 
    plot_grid(top_grid, NULL, bottom_grid, ncol = 1, rel_heights = c(.25, .025, 1)) +
    theme(panel.background = element_rect(color = "white", fill = "white"))

ggsave("./pdf/fig4.pdf", 
       final_grid,
       width = 179,
       height = 190,
       units = "mm",
       dpi = 600,
       device = cairo_pdf
       )
