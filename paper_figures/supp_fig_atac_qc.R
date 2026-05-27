library(DESeq2)
library(tidyverse)
library(glue)

# meta data
samplesheet <- 
    read_csv("../atacseq/samplesheet.csv") |>
    mutate(f1 = basename(fastq_1),
           donor_id = str_extract(f1, "20221025_([^_]+).*", group = 1),
           replic = paste0("REP", replicate)) |>
    select(condition = sample, replic, donor_id) |>
    distinct()

col_data <-
    samplesheet |>
    unite("sample_id", c(condition, replic), sep = "_", remove = FALSE) |>
    mutate(sample_type = case_when(donor_id == "3donors" ~ "pool",
                                   .default = "donor")) |>
    group_by(condition) |>
    filter(any(sample_type == "pool") & !all(sample_type == "pool")) |>
    ungroup() |>
    mutate(condition = recode(condition, 
                              "unst_0" = "Unstim 0h", 
                              "IL4_24" = "IL-4c 24h",
                              "TLR7_24" = "TLR7c 24h",
                              "BCR_24" = "BCRc 24h",
                              "DN2_24" = "DN2c 24h"),
           condition = factor(condition, levels = c("Unstim 0h", "IL-4c 24h", "TLR7c 24h", "BCRc 24h", "DN2c 24h"))) |>
    arrange(condition, replic) |>
    column_to_rownames("sample_id")

# Count data
dat <-
    "../atacseq/results/bwa/merged_replicate/macs2/narrow_peak/consensus/consensus_peaks.mRp.clN.featureCounts.txt" |>
    read_table(skip = 1)

colnames(dat) <- gsub("\\.mLb\\.\\clN\\.sorted\\.bam", "", colnames(dat))

# Interval info
interval_table <-
    dat |>
    select(Geneid:Length) |>
    column_to_rownames("Geneid")

#counts
count_table <- 
    dat |>
    select(Geneid, all_of(rownames(col_data))) |>
    column_to_rownames("Geneid")

## RUN DESEQ2
dds <- 
    DESeqDataSetFromMatrix(countData = count_table, 
                           colData = col_data, 
                           design = ~ 1)

dds <- estimateSizeFactors(dds)
keep_intervals <- rowSums(counts(dds) >= 10) >= 4
dds_filtered <- dds[keep_intervals, ]

norm_counts <- counts(dds_filtered, normalized = TRUE)

rlog_out <- rlog(dds_filtered, blind = TRUE)
rlog_mat <- assay(rlog_out)

rlog_df <- 
    as.data.frame(rlog_mat) |>
    as_tibble(rownames = "interval_id") |>
    pivot_longer(-interval_id, names_to = "sample_id") |>
    left_join(as_tibble(col_data, rownames = "sample_id"), join_by(sample_id)) |>
    select(interval_id, condition, replic, value) |>
    pivot_wider(names_from = replic, values_from = value) |>
    pivot_longer(REP2:REP4, names_to = "donor")

corrs_df <- 
    rlog_df |>
    group_by(condition, donor) |>
    summarize(pearson_r = cor(value, REP1, method = "pearson")) |>
    ungroup() |>
    mutate(label = paste0("R = ", round(pearson_r, 2)))

plot_out <- 
    ggplot(rlog_df, aes(x = REP1, y = value)) +
    geom_bin2d(bins = 100, show.legend = FALSE) + 
    scale_fill_viridis_c() +
    geom_text(
        data = corrs_df,
        aes(label = label, x = -Inf, y = Inf),
        hjust = -0.2, 
        vjust = 1.5, 
        inherit.aes = FALSE, 
        size = 3, fontface = "bold"
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_grid(donor ~ condition) +
    labs(
        x = "Pool Expression (rlog)",
        y = "Individual Donor Expression (rlog)"
    ) +
    theme_bw() +
    theme(
        legend.text  = element_text(size = 7),
        legend.title = element_text(size = 8),
        strip.background = element_rect(fill = "grey90")
    )

# plot QC metrics
stim_colors <- 
    read_tsv("./figure_colors.txt", col_names = c("stim", "timep", "col")) |>
    mutate(condition = glue("{stim} {timep}hrs")) |>
    select(condition, col) |>
    deframe()

frip_df <- 
    "../atacseq/results/multiqc/narrow_peak/B_cell_atac_nfcore_multiqc_report_data/multiqc_mlib_frip_score-plot.txt" |>
    read_tsv() |>
    pivot_longer(-Sample) |>
    drop_na() |>
    select(Sample, frip = value) |>
    left_join(unite(samplesheet, "Sample", c(condition, replic), sep = "_", remove = FALSE), join_by(Sample)) |>
    filter(donor_id != "3donors") |>
    select(sample_id = Sample, condition, donor_id, frip) |>
    separate(condition, c("stim", "timep"), sep = "_") |>
    mutate(stim = recode(stim, 
			 "unst" = "Unstim", 
			 "IL4" = "IL-4c",
			 "TLR7" = "TLR7c",
			 "BCR" = "BCRc",
			 "DN2" = "DN2c"),
           stim = factor(stim, levels = c("Unstim", "IL-4c", "TLR7c", "BCRc", "DN2c"))) |> 
    arrange(stim, as.numeric(timep), donor_id) |>
    mutate(condition = glue("{stim} {timep}hrs"),
	   condition = fct_inorder(condition))

frip_plot <- 
    ggplot(data = frip_df,
	   aes(x = frip, y = donor_id)) +
    geom_col(alpha = .75, color = "black", linewidth = .2) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap(~condition, ncol = 1, strip.position = "right") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 7),
	  axis.text.y = element_blank(),
	  axis.title = element_text(size = 8),
	  strip.text.y.right = element_text(size = 8, angle = 0, hjust = 0),
	  panel.spacing = unit(0.1, "lines"),
	  panel.grid.minor.x = element_blank(),
	  plot.title = element_text(size = 8),
	  plot.background = element_rect(fill = "white", color = "white")
    ) +
    labs(y = "Donor", x = "FriP", title = "Fraction of reads in peaks") 

peak_annot <- 
    "../atacseq/results/multiqc/narrow_peak/B_cell_atac_nfcore_multiqc_report_data/multiqc_mlib_peak_annotation-plot.txt" |>
    read_tsv() |>
    pivot_longer(-Sample, names_to = "annotation") |>
    left_join(unite(samplesheet, "Sample", c(condition, replic), sep = "_", remove = FALSE), join_by(Sample)) |>
    filter(donor_id != "3donors") |>
    separate(condition, c("stim", "timep"), sep = "_") |>
    mutate(stim = recode(stim, 
			 "unst" = "Unstim", 
			 "IL4" = "IL-4c",
			 "TLR7" = "TLR7c",
			 "BCR" = "BCRc",
			 "DN2" = "DN2c"),
           stim = factor(stim, levels = c("Unstim", "IL-4c", "TLR7c", "BCRc", "DN2c"))) |> 
    arrange(stim, as.numeric(timep), donor_id) |>
    mutate(condition = glue("{stim} {timep}hrs"),
	   condition = fct_inorder(condition)) |>
    select(donor_id, condition, annotation, value)
   
annot_colors <- 
    c(
      "promoter-TSS" = "#E41A1C",  # Red
      "exon"         = "#377EB8",  # Blue
      "intron"       = "#A6CEE3",  # Light Blue
      "TTS"          = "#FF7F00",  # Orange
      "Intergenic"   = "#999999",  # Medium Grey
      "Unassigned"   = "#525252"   # Dark Grey
    )

annot_plot <- 
    ggplot(data = peak_annot) +
    geom_col(aes(x = value, y = donor_id, fill = annotation), position = "fill") +
    scale_x_continuous(labels = scales::percent,
		       expand = c(0, 0)) +
    scale_fill_manual(values = annot_colors) + 
    facet_wrap(~condition, ncol = 1, strip.position = "right") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 7),
	  axis.text.y = element_blank(),
	  axis.title = element_text(size = 8),
	  strip.text.y.right = element_text(size = 8, angle = 0, hjust = 0),
	  panel.spacing = unit(0.1, "lines"),
	  panel.grid.minor.x = element_blank(),
	  legend.position = "right",
	  legend.text  = element_text(size = 7),
	  legend.title = element_text(size = 8),
	  legend.key.size = unit(0.3, "cm"),
	  legend.box.margin = margin(0, 0, 0, -10), 
	  legend.spacing.y = unit(0.2, "cm"),
	  legend.key = element_blank(),
	  plot.title = element_text(size = 8),
	  plot.background = element_rect(fill = "white", color = "white")
    ) +
    labs(y = "Donor", x = "Percentage", title = "HOMER peak annotation") 

library(patchwork)

ggsave("atac_qc.png", 
       free(plot_out) / (frip_plot + annot_plot) +
           plot_layout(heights = c(2, 1)) +
           plot_annotation(tag_levels = "A"),
       height = 6.5, width = 6.5)



