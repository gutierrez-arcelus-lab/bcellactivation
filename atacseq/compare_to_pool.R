library(DESeq2)
library(tidyverse)
library(glue)

# meta data
samplesheet <- 
    read_csv("./samplesheet.csv") |>
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
    "./results/bwa/merged_replicate/macs2/narrow_peak/consensus/consensus_peaks.mRp.clN.featureCounts.txt" |>
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
        title = "Replicate Consistency: Individual Donors vs. Independent Pool",
        x = "Pool Expression (rlog)",
        y = "Individual Donor Expression (rlog)"
    ) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "grey90"))

ggsave("./donors_vs_pool.png", plot_out, width = 6.5, height = 4, dpi = 300)
