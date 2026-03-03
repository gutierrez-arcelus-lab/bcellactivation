library(tidyverse)
library(ggrepel)
library(scico)
library(cowplot)
library(glue)
library(tidytext)

# colors
stim_colors <- 
    read_tsv("./figure_colors.txt", col_names = c("stim", "time", "color")) |>
    unite("stim", c(stim, time), sep = " ") |>
    filter(grepl("0$|24$", stim), 
	   grepl("^Unstim\\s|^IL-4c\\s|^TLR7c\\s|^BCRc\\s|^DN2c\\s", stim)) |>
    mutate(stim = paste0(stim, "h")) |>
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
	   stim = str_replace(stim, "IL4", "IL-4c"),
	   stim = str_replace(stim, "BCR", "BCRc"),
	   stim = str_replace(stim, "TLR7", "TLR7c"),
	   stim = str_replace(stim, "DN2", "DN2c"),
	   stim = paste0(stim, "h")) |>
    select(stim, donor_id, donor = replicate) |>
    distinct() |>
    filter(donor_id != "3donors")

# Figure A ####################################################################
pca_file <- "../atacseq/results_deseq2/pcadata_5000peaks.rds"
#pca_file <- "../atacseq/results_deseq2/mRp/pcadata_5000peaks.rds"

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
    geom_point(aes(fill = stim), shape = 21, size = 2.5, stroke = .5) +
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
	  legend.key.height = unit(.1, "in"), 
	  legend.box.spacing = unit(1, "pt")
	  ) +
    coord_cartesian(clip = "off") +
    labs(x = sprintf("PC1 (%s%%)", perc_var[1]),
	 y = sprintf("PC2 (%s%%)", perc_var[2]),
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

fig_a_grid <- plot_grid(fig_a_title, 
			plot_grid(NULL, fig_a, NULL, nrow = 1, rel_widths = c(0.05, 1, 0.05)), 
			ncol = 1, rel_heights = c(.1, 1))


# Figure B ####################################################################
#da_files <- 
#    list.files("../atacseq/results_deseq2/mRp",
#	       pattern = ".+vs.+\\.tsv$",
#	       full.names = TRUE)

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
    theme_minimal() + 
    theme(
	  panel.grid = element_blank(),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.ticks.y = element_blank(),
	  legend.position = "none",
	  plot.margin = unit(c(0, 0, 0, 0), "lines")
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
    coord_cartesian(xlim = c(0.5, 1.5), expand = FALSE)

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
    coord_cartesian(ylim = c(0.5, 1.5), expand = FALSE)

fig_b <- 
    plot_grid(NULL, cond_horizontal, cond_vertical, da_plot,
	      ncol = 2, nrow = 2, align = "hv",
	      rel_heights = c(.125, 1), rel_widths = c(0.05, 1))

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
    plot_grid(fig_b_title, NULL, fig_b, ncol = 1, rel_heights = c(.1, .01, 1)) +
    theme(plot.margin = margin(r = 5.5 * 2, unit = "pt"))
	

# Figure C ####################################################################
 
## Filter motifs for TFs that are expressed in the RNA-seq data
motif_gene_mappings <- 
    "./data/homer_motif_gene_mappings.xlsx" |>
    readxl::read_excel() |>
    select(motif_info = 1, gene_name = 2, notes = 3) |>
    mutate(gene_name = case_match(gene_name,
				  "STAT5A, STAT5B" ~ "STAT5A/STAT5B",
				  "CREB1, ATF1" ~ "CREB1/ATF1",
				  "IRF4, IRF8, BATF" ~ "IRF8, BATF",
				  "NFATC1, NFATC2" ~ "NFATC1/NFATC2",
				  "FOS, JUN, FOSL1/2" ~ "FOS, JUN, FOSL1/FOSL2",
				  .default = gene_name))

gene_expression <-
    "../bcell_lowinput/data/expression_pooled_reps.tsv" |>
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
    group_by(gene_id, gene_name, stim, timep) |>
    summarize(mean_tpm = mean(tpm)) |>
    ungroup() |>
    group_by(gene_id, gene_name) |>
    filter(any(mean_tpm >= 5)) |>
    ungroup() |>
    distinct(gene_id, gene_name) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

stims <- read_lines("../atacseq/homer/stims.txt")

homer_df <-
    glue("../atacseq/homer/results_vert/{stims}/knownResults.txt") |>
    #glue("../atacseq/homer/results_vert_mRp/{stims}/knownResults.txt") |>
    setNames(paste0(stims, "c")) |>
    map_dfr(~read_tsv(.) |>
            janitor::clean_names() |>
            {function(x) setNames(x, str_remove(names(x), "_of_\\d+$"))}(),
            .id = "stim") |>
    rename(motif_info = motif_name) |>
    mutate_at(vars(starts_with("percent_of_")), parse_number) |>
    mutate(stim = factor(stim, levels = paste0(stims, "c")),
           fc = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) |>
    extract(motif_info,
            c("motif", "motif_family", "dataset", "db"),
            "(.+?)(?:\\((.+)\\))?/([^/]+)/(.+)",
            remove = FALSE) |>
    mutate(log10p = log_p_value/log(10)) |>
    select(stim, motif_info, motif, motif_family, dataset, consensus,
           pct_target =  percent_of_target_sequences_with_motif,
           pct_bg = percent_of_background_sequences_with_motif,
           fc, log10p, q_value = q_value_benjamini)

get_expression_status <- function(g) {

    if ( !grepl(",|/", g) ) {    
	
	res <- g %in% genes_expressed$gene_name

    } else if ( grepl("/", g) && !grepl(",", g) ) {
	
	res <- any(unlist(str_split(g, "/")) %in% genes_expressed$gene_name)

    } else if ( grepl(",", g) && !grepl("/", g) ) { 
       
	res <- all(unlist(str_split(g, ", ")) %in% genes_expressed$gene_name)

    } else if ( grepl("/", g) && grepl(",", g) ) {

	g1 <- unlist(str_split(g , ", ")) |> {function(x) keep(x, !grepl("/", x))}()
	g2 <- unlist(str_split(g , ", ")) |> {function(x) keep(x, grepl("/", x))}() |> str_split("/") |> unlist()  

	res <- all(g1 %in% genes_expressed$gene_name) && any(g2 %in% genes_expressed$gene_name)
    }
}

motifs_expressed_tf <- 
    homer_df |>
    filter(q_value <= 0.01) |>
    distinct(motif_info) |>
    left_join(motif_gene_mappings) |>
    mutate(is_expressed = map_lgl(gene_name, get_expression_status)) |>
    filter(is_expressed) |>
    select(motif_info)

homer_top <- 
    homer_df |>
    inner_join(motifs_expressed_tf) |>
    group_by(motif) |>
    slice_min(log10p) |>
    ungroup() |>
    top_n(50, -log10p)

motif_family_order <- 
    homer_top |>
    group_by(motif_family) |>
    slice_min(log10p) |>
    arrange(log10p) |>
    select(motif_info, motif_family, fc, log10p, q_value) |>
    pull(motif_family)

homer_top_inall <- 
    homer_df |>
    inner_join(distinct(homer_top, motif, motif_family)) |>
    arrange(log10p) |>
    mutate(motif = fct_inorder(motif),
	   motif = fct_rev(motif),
	   motif_family = factor(motif_family, levels = motif_family_order))

fig_c <- 
    ggplot(homer_top_inall, 
	   aes(x = stim, y = motif)) +
    geom_point(aes(fill = log2(fc), size = -log10p), 
	       shape = 21, stroke = .2) +
    geom_text(data = filter(homer_top_inall, q_value <= 0.01), 
	      aes(x = stim, y = motif, label = "*"), 
	      size = 8, fontface = "bold", size.unit = "pt", 
	      nudge_x = .35, nudge_y = -.2) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "left") +
    scale_fill_gradient2(high = '#e66101', mid = 'white', low = '#5e3c99', midpoint = 0) +
    scale_size(range = c(2, 5)) +
    facet_grid(rows = vars(motif_family), scales = "free", space = "free", switch = "x") +
    theme_minimal() +
    theme(
	  axis.text.x.top = element_text(size = 8, hjust = 0.5), 
	  axis.text.y = element_text(size = 8), 
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


# Figure D ####################################################################
traits <- 
    read_tsv("../atacseq/ldsc/data/traits.txt", col_names = c("directory", "trait", "gwas", "ref")) |>
    select(gwas, trait)

results <- 
    read_tsv("../atacseq/ldsc/compiled_results.tsv") |>
    #read_tsv("../atacseq/ldsc/mRp/compiled_results.tsv") |>
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
    geom_text(data = filter(results, pfdr <= 1.01),
	      aes(x = set, y = gwas, label = "*"), 
	      size = 10, fontface = "bold", size.unit = "pt", nudge_y = -.2) +
    scale_x_discrete(position = "top") +
    scale_fill_gradient2(high = '#e66101', mid = 'white', low = '#5e3c99', midpoint = 0) +
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
	      labels = c("a", "b"), label_size = 12)

bottom_grid <- 
    plot_grid(fig_c_grid, NULL, fig_d_grid, nrow = 1, rel_widths = c(.85, .0, 1),
	      labels = c("c", "", "d"), label_size = 12)

final_grid <- 
    plot_grid(top_grid, NULL, bottom_grid, ncol = 1, rel_heights = c(.25, .01, 1)) +
    theme(panel.background = element_rect(color = "white", fill = "white"))

ggsave("fig4.png", final_grid, width = 6.5, height = 7.5, dpi = 300)
ggsave("./high_res/fig4.png", final_grid, width = 6.5, height = 7.5, dpi = 600)

