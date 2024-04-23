library(tidyverse)
library(tidytext)
library(broom)
library(cowplot)
library(extrafont)

stim_colors <- 
    c("unstday0" = "grey",
      "BCR" = "#4DBBD5FF",
      "TLR7" = "#91D1C2FF",
      "DN2" = "#E64B35FF")

stim_shades <- 
    c("unstday0_raw_reads" = "black",
      "unstday0_Unique" = "grey60",
      "BCR_raw_reads" = "midnightblue",
      "BCR_Unique" = "#94a4db",
      "TLR7_raw_reads" = "forestgreen",
      "TLR7_Unique" = "#8dbf7d",
      "DN2_raw_reads" = "#801c05",
      "DN2_Unique" = "#d08d74")

# QC
temp_dir <- system("echo $TEMP_WORK", intern = TRUE)

read_fastq_counts <- function(qc_dir) {

    file.path(qc_dir, "mqc_fastqc_sequence_counts_plot_1.txt") |>
    read_tsv() |>
    extract(Sample, 
	    c("sample_id", "stim", "r"), 
	    "\\d+_(\\d+[_123]*)_(.+)_MG\\d+_.+_(R[12])") |>
    mutate(donor_id = str_extract(sample_id, "[^_]+")) |>
    pivot_longer(`Unique Reads`:`Duplicate Reads`, names_to = "read_type", values_to = "n") |>
    mutate(read_type = sub("\\sReads$", "", read_type)) |>
    select(donor_id, sample_id, stim, r, read_type, n)
}

qc_df <- 
    file.path(temp_dir, 
	      "fastq/highinput/multiqc_data") |>
    read_fastq_counts()

# Total raw reads
raw_df <- 
    file.path(temp_dir, "fastq/highinput_raw/multiqc_data/multiqc_general_stats.txt") |>
    read_tsv() |>
    select(Sample, raw_reads = ends_with("total_sequences")) |>
    extract(Sample, 
	    c("sample_id", "stim", "r"), 
	    "\\d+_(\\d+[_123]*)_(.+)_MG\\d+_.+_(R[12])") |>
    mutate(donor_id = str_extract(sample_id, "[^_]+")) |>
    filter(r == "R1") |>
    select(donor_id, sample_id, stim, raw_reads)

qc_df_filter <- qc_df |>
    filter(r == "R1") |>
    select(-r) 

plot_df <- raw_df |>
    pivot_longer(raw_reads, names_to = "read_type", values_to = "n") |>
    bind_rows(qc_df_filter) |>
    mutate(stim = recode(stim, "TR7" = "TLR7"),
	   stim = factor(stim, levels = names(stim_colors))) |>
    arrange(sample_id, stim, read_type) |>
    filter(read_type != "Duplicate") |>
    mutate(lab = paste(stim, read_type, sep = "_"),
	   read_type = factor(read_type, levels = c("raw_reads", "Unique")))

total_reads <- 
    ggplot(plot_df, 
	   aes(x = n,
	       y = reorder_within(sample_id, by = n, within = stim, fun = "max"))) +
    geom_col(aes(fill = lab),
	     position = "dodge",
	     color = "black", linewidth = .1, width = 1) +
    scale_fill_manual(values = stim_shades) +
    scale_x_continuous(labels = function(x) x/1e6L, 
		       breaks = seq(0, 70e6, 10e6),
		       expand = c(0, 0)) +
    scale_y_reordered() +
    facet_wrap(~stim, scales = "free", ncol = 2) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "grey90", size = .5),
          panel.grid.major.y = element_blank(),
	  panel.spacing = unit(2, "lines"),
	  axis.title = element_text(size = 16),
	  strip.text = element_text(size = 16),
	  axis.text.x = element_text(size = 12),
	  axis.text.y = element_text(hjust = 0),
          plot.background = element_rect(fill = "white", color = "white"),
          legend.position = "none") +
    labs(x = "Million reads", y = NULL, fill = "Stim")

ggsave("./plots/total_reads.png", total_reads, width = 8, height = 8)


stim_colors2 <- 
    c("unstday0" = "#898e9f",
      "BCR" = "#003967",
      "TLR7" = "#637b31",
      "DN2" = "#a82203")

total_reads_raw_plot <- 
    raw_df |>
    mutate(stim = recode(stim, "TR7" = "TLR7"),
	   stim = factor(stim, levels = names(stim_colors2))) |>
    ggplot(aes(x = raw_reads,
	       y = reorder_within(sample_id, by = raw_reads, within = stim))) +
    geom_col(aes(fill = stim),
	     color = "black", linewidth = .1, width = 1) +
    scale_fill_manual(values = stim_colors2) +
    scale_x_continuous(labels = function(x) x/1e6L,
		       limits = c(0, max(raw_df$raw_reads)),
		       expand = c(0, 0)) +
    scale_y_reordered() +
    facet_wrap(~stim, scales = "free", ncol = 2) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
	  panel.grid.major.x = element_line(color = "grey90", linewidth = .5),
	  panel.grid.major.y = element_blank(),
	  panel.spacing = unit(2, "lines"),
	  axis.title = element_text(size = 16, family = "Arial"),
	  axis.text.x = element_text(size = 12, family = "Arial"),
	  axis.text.y = element_text(hjust = 0, family = "Arial"),
	  strip.text = element_text(size = 16, family = "Arial"),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = "Total reads (Million)", y = NULL) +
    guides(fill = "none")

ggsave("./plots/total_reads_raw.png", total_reads_raw_plot, width = 8, height = 8)

raw_df |>
    select(sample_id, stim, n = raw_reads) |>
    filter(n < 40e6L) |>
    mutate(stim = factor(stim, levels = names(stim_colors2))) |>
    arrange(stim, n) |>
    unite("sample_id", c(sample_id, stim), sep = "_") |>
    write_tsv("low_total_reads_samples_mga.tsv")

# Resequencing of low depth samples
reseq_samples <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/resequencing/231227_MG10430_fastq/" |>
    list.files(pattern = "fastq\\.gz") |>
    {function(x) tibble(f = x)}() |>
    extract(f, c("sample_id", "stim", "r"), "\\d+_(\\d+[_123]*)_(.+)_MG\\d+_.+_(R[12])") |>
    mutate(i = 1) |>
    pivot_wider(names_from = r, values_from = i) |>
    unite("sample_id", c(sample_id, stim), sep = "_")

reseq_samples |>
    filter(is.na(R1) | is.na(R2)) |>
    mutate_at(vars(R1, R2), ~replace_na(., 0))

raw_rsq_df |>
    select(sample_id, stim) |>
    unite("sample_id", c(sample_id, stim), sep = "_") |>
    anti_join(reseq_samples, join_by(sample_id))



# Additional resequencing
raw_rsq_df <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq",
	      "resequencing/240122_MG10430_fastq/multiqc_data",
	      "multiqc_general_stats.txt") |>
    read_tsv() |>
    select(Sample, raw_reads = ends_with("total_sequences")) |>
    extract(Sample, 
	    c("sample_id", "stim", "r"), 
	    "\\d+_(\\d+[-123]*)-(.+)_MG\\d+_.+_(R[12])") |>
    mutate(sample_id = sub("-", "_", sample_id),
	   donor_id = str_extract(sample_id, "[^_]+")) |>
    filter(r == "R1") |>
    select(donor_id, sample_id, stim, raw_reads)

raw_all_df <-
    bind_rows(raw_df, raw_rsq_df) |>
    group_by(sample_id, stim) |>
    summarise_at(vars(raw_reads), sum) |>
    ungroup() |>
    mutate(stim = recode(stim, "TR7" = "TLR7"),
	   stim = factor(stim, levels = names(stim_colors2)))

total_reads_raw_rsq_plot <- 
    raw_all_df |>
    ggplot(aes(x = raw_reads,
	       y = reorder_within(sample_id, by = raw_reads, within = stim))) +
    geom_col(aes(fill = stim),
	     color = "black", linewidth = .1, width = 1) +
    scale_fill_manual(values = stim_colors2) +
    scale_x_continuous(labels = function(x) x/1e6L,
		       limits = c(0, max(raw_df$raw_reads)),
		       expand = c(0, 0)) +
    scale_y_reordered() +
    facet_wrap(~stim, scales = "free", ncol = 2) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
	  panel.grid.major.x = element_line(color = "grey90", linewidth = .5),
	  panel.grid.major.y = element_blank(),
	  panel.spacing = unit(2, "lines"),
	  axis.title = element_text(size = 16, family = "Arial"),
	  axis.text.x = element_text(size = 12, family = "Arial"),
	  axis.text.y = element_text(hjust = 0, family = "Arial"),
	  strip.text = element_text(size = 16, family = "Arial"),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = "Total reads (Million)", y = NULL) +
    guides(fill = "none")

ggsave("./plots/total_reads_raw_rsq.png", total_reads_raw_rsq_plot, width = 8, height = 8)





# RNA extraction data
library(readxl)

rna_extraction_df <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/DNA-RNAextraction_MGB.xlsx" |>
    read_excel(skip = 2) |>
    janitor::clean_names() |>
    drop_na(date) |>
    select(sample_id, n_cells, rna_conc = 4, rna_total = 5, dna_conc = 6, dna_total = 7) |>
    mutate(sample_id = sub("TR7", "TLR7", sample_id),
	   n_cells = sub("million", "e+06", n_cells),
	   n_cells = sub(",", ".", n_cells),
	   n_cells = sub(" ", "", n_cells),
	   n_cells = as.integer(n_cells)) |>
    mutate_at(vars(rna_conc:dna_total), as.numeric) |>
    extract(sample_id, c("sample_id", "stim"), "(\\d+[_123]*)_(.+)") |>
    mutate(stim = factor(stim, levels = names(stim_colors)))


plot_df2 <- 
    left_join(plot_df, rna_extraction_df, join_by(sample_id, stim)) |>
    pivot_wider(names_from = read_type, values_from = n) |>
    mutate(Total = Unique + Duplicate) |>
    select(-Duplicate) |>
    pivot_longer(c(Unique, Total), names_to = "read_type", values_to = "n_reads")
    

scatter_p <- ggplot(plot_df2, aes(rna_conc, y = n_reads)) +
    geom_point(aes(color = stim), size = 2) +
    scale_y_continuous(labels = function(x) x/1e6L) +
    scale_color_manual(values = stim_colors) +
    facet_wrap(~read_type, ncol = 1, scales = "free") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank()) +
    labs(x = "RNA concentration (ng/ul)", y = "Million reads", color = "Stim:")
    
ggsave("./plots/reads_by_rna.png", scatter_p, height = 5, width = 5)



# PCA 
quant <- read_tsv("./results/salmon_genes.tsv")

variable_genes <- quant |>
    group_by(gene_id, gene_name) |>
    summarise(v = var(tpm)) |>
    ungroup() |>
    top_n(5000, v) |>
    arrange(desc(v))

variable_matrix <- quant |>
    filter(gene_id %in% variable_genes$gene_id) |>
    select(-gene_name, -counts) |>
    pivot_wider(names_from = gene_id, values_from = tpm) |>
    unite("id", c(sample_id, stim), sep = "_") |>
    column_to_rownames("id")

pca <- prcomp(variable_matrix, center = TRUE, scale. = TRUE, rank. = 10)

var_exp <- pca |>
    tidy("pcs") |>
    filter(PC %in% 1:10) |>
    ggplot(aes(x = factor(PC), y = percent)) +
	geom_col(alpha = 0.7) +
	scale_y_continuous(labels = scales::label_percent()) +
    theme_minimal() +
    theme(text = element_text(size = 14),
	  panel.grid.major.x = element_blank()) +
    labs(x = "PC", y = "Variance explained")

pc_scores <- as_tibble(pca$x, rownames = "id")

pca_df <- pc_scores |>
    separate(id, c("sample_id", "stim"), sep = "_") |>
    mutate(stim = factor(stim, levels = names(stim_colors)))

pca_plot <- ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(fill = stim), alpha = .25, size = 4, shape = 21) +
    ggrepel::geom_text_repel(aes(label = sample_id, color = stim),
			     fontface = "bold", size = 2.5, 
			     max.overlaps = 1000, min.segment.length = 0,
			     show.legend = FALSE) +
    scale_color_manual(values = stim_colors) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(text = element_text(size = 14),
	  panel.grid = element_blank(),
          legend.key.height = unit(.25, "cm"),
	  legend.position = c(.9, .1)) +
    guides(fill = guide_legend(ncol = 1))

pca_out <- plot_grid(var_exp, pca_plot, nrow = 1, rel_widths = c(.5, 1)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/pca.png", pca_out, width = 8, height = 4)

ggsave("./plots/pca_labels.png", pca_plot + theme(plot.background = element_rect(fill = "white", color = "white")))


# Edger
cpm_df <- read_tsv("./edger_cpm.tsv") |>
    mutate(stim_test = factor(stim_test, levels = c("BCR", "TLR7", "DN2")),
	   stim = factor(stim, levels = c("unstday0", "BCR", "TLR7", "DN2")))

edger_res <- read_tsv("./edger_results.tsv") |>
    mutate(stim_test = factor(stim_test, levels = c("BCR", "TLR7", "DN2")))

fold_change_summ <- 
    edger_res |> 
    group_by(stim = stim_test, is_up = logFC >= 0) |>
    summarise(n = sum(abs(logFC) >= 1 & FDR < .05)) |>
    ungroup() |>
    mutate(is_up = recode(as.character(is_up), "FALSE" = "Down", "TRUE" = "Up")) |>
    ggplot(aes(x = is_up, y = n, fill = stim)) +
	geom_col() +
	scale_fill_manual(values = c("BCR" = "#6996e3", "TLR7" = "#748f46", "DN2" = "#d76b51")) +
	facet_wrap(~stim, nrow = 1) +
	theme_minimal() +
	theme(legend.position = "none",
	      panel.grid.major.y = element_blank(),
	      plot.background = element_rect(color = "white", fill = "white")) +
	labs(x = "Up/Down regulated", y = "Number of genes")

ggsave("./plots/fold_change_summary.png", 
       fold_change_summ, 
       width = 4, height = 3, dpi = 600)



sle_genes <- read_tsv("../bcell_scrna/reported_genes.tsv")

plot_fc <- 
    edger_res |> 
    filter(gene_name %in% sle_genes$gene, FDR < 0.05) |>
    group_by(stim_test) |>
    top_n(25, abs(logFC)) |>
    ungroup() |>
    ggplot(aes(x = logFC,y = reorder_within(gene_name, within = stim_test, by = logFC))) +
	geom_col(aes(fill = stim_test)) +
	scale_y_reordered() +
	scale_fill_manual(values = c("BCR" = "#6996e3", "TLR7" = "#748f46", "DN2" = "#d76b51")) +
	facet_wrap(~stim_test, nrow = 1, scales = "free") +
	theme_bw() +
	theme(text = element_text(size = 12),
	      legend.position = "none",
	      panel.grid.minor.x = element_blank()) +
	labs(x = "log2 Fold-change", y = NULL)

ggsave("./plots/fold_change.png", plot_fc, width = 10, height = 5)


stat1 <- cpm_df |> 
    filter(gene_name == "STAT1") |>
    ggplot(aes(x = stim, y = cpm)) +
	geom_line(aes(group = sample_id), linewidth = .2) +
	geom_point(aes(fill = stim), size = 4, shape = 21, stroke = .2) +
	scale_fill_manual(values = stim_colors) +
	facet_wrap(~stim_test, scales = "free_x", nrow = 1) +
	theme_bw() +
	theme(legend.position = "none",
	      panel.grid.minor.y = element_blank()) +
	labs(x = NULL, y = "Counts per Million")

ggsave("./plots/stat1.png", stat1, height = 3)


irf8 <- cpm_df |> 
    filter(gene_name == "IRF8") |>
    ggplot(aes(x = stim, y = cpm)) +
	geom_line(aes(group = sample_id), linewidth = .2) +
	geom_point(aes(fill = stim), size = 4, shape = 21, stroke = .2) +
	scale_fill_manual(values = stim_colors) +
	facet_wrap(~stim_test, scales = "free_x", nrow = 1) +
	theme_bw() +
	theme(legend.position = "none",
	      panel.grid.minor.y = element_blank()) +
	labs(x = NULL, y = "Counts per Million")

ggsave("./plots/irf8.png", irf8, height = 3)

# GO

go_res <- read_tsv("./top_diff_go.tsv") |>
    mutate(stim = factor(stim, levels = names(stim_colors[-1]))) |>
    arrange(desc(-log10(P.Up))) |>
    mutate(Term = fct_inorder(Term))

goplot <- ggplot(go_res) +
    geom_hline(yintercept = -log10(0.05), linewidth = .5) +
    geom_hline(yintercept = -log10(0.01), linewidth = .5, linetype = "dashed") +
    geom_col(aes(x = stim, y = -log10(P.Up), fill = stim), position = "dodge") +
    scale_fill_manual(values = stim_colors[-1]) +
    facet_wrap(~Term, nrow = 4,
	       labeller = labeller(Term = label_wrap_gen(width = 25))) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_shape_manual(values = 1:nlevels(mbv_df$donor_vcf) - 1) +
    labs(x = NULL)

ggsave("./plots/go.png", goplot, width = 10, height = 6)


# MBV analysis

mbv_files <- list.files("./results/mbv/")

mbv_df <- file.path("./results/mbv", mbv_files) |>
    setNames(mbv_files) |>
    map_dfr(~read_delim(., delim = " "), .id = "bam") |>
    extract(bam, c("donor_id_bam", "bam_stim"), "(\\d+\\.\\d)_([^.]+)\\.txt") |>
    separate(donor_id_bam, c("donor_bam", "rep_bam"), sep = "\\.", remove = FALSE) |>
    extract(SampleID, c("donor_vcf"), "\\d+_[^-]+-(\\d+)")

donor_id_order <- mbv_df |> 
    distinct(donor_id_bam, bam_stim) |> 
    count(donor_id_bam) |>
    arrange(desc(n), donor_id_bam) |>
    pull(donor_id_bam)

mbv_df <- 
    mbv_df |>
    mutate(donor_id_bam = factor(donor_id_bam, levels = donor_id_order),
	   bam_stim = factor(bam_stim, levels = c("unstday0", "TLR7", "BCR", "DN2")), 
	   id_match = donor_bam == donor_vcf,
	   lab = ifelse(id_match == FALSE & perc_hom_consistent > .8 & perc_het_consistent > .8,
			donor_vcf,
			NA))

mbv_plot <- 
    ggplot(mbv_df, 
       aes(perc_het_consistent, perc_hom_consistent)) +
    geom_point(aes(color = id_match), size = 2, shape = 1, stroke = 2) +
    ggrepel::geom_text_repel(aes(label = lab), 
			     size = 3,
			     min.segment.length = 0,
			     fontface = "bold") +
    scale_color_manual(values = c("TRUE" = "midnightblue", "FALSE" = "tomato2")) +
    scale_x_continuous(limits = c(0, 1), 
		       breaks = seq(0, 1, .5), 
		       labels = c("0", "0.5", "1")) + 
    scale_y_continuous(limits = c(0, 1), 
		       breaks = seq(0, 1, .5),
		       labels = c("0", "0.5", "1")) +
    facet_wrap(donor_id_bam~bam_stim, ncol = 8) +
    theme_bw() +
    theme(text = element_text(size = 12),
	  panel.grid = element_blank(),
	  legend.position = c(.9, .025)) +
    labs(color = "Sample ID\nmatch",
	 x = "Fraction concordant heterozygous sites",
	 y = "Fraction concordant homozygous sites")

ggsave("./plots/mbv.png", mbv_plot, width = 10, height = 12)

mbv_df |>
    group_by(donor_id_bam) |>
    filter(!any(perc_hom_consistent > .8 & perc_het_consistent > .8)) |>
    select(donor_id_bam, bam_stim:id_match) |>
    print(n = Inf, width = Inf)

mbv_df |>
    filter(donor_bam == "10051708") |>
    print(width = Inf)
