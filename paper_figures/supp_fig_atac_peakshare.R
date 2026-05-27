library(tidyverse)
library(patchwork)
library(glue)

stims <- c("IL4", "TLR7", "BCR", "DN2")

stim_colors <- 
    read_tsv("./figure_colors.txt", col_names = c("stim", "timep", "col")) |>
    filter(timep == 72) |>
    filter(stim %in% c("IL-4c", "TLR7c", "BCRc", "DN2c"))

da_data <- 
    glue("../atacseq/results_deseq2/{stims}_24vsunst_24.tsv") |>
    setNames(stims) |>
    map_dfr(read_tsv, .id = "stim") |>
    mutate(stim = factor(stim, levels = stims))

peak_summary <- 
    da_data |>
    group_by(stim) |>
    summarise(Total = n(),
              Tested = sum(!is.na(padj))) |>
    ungroup() |>
    pivot_longer(-stim)

peak_lists <- 
    da_data |>
    filter(padj <= 0.01, log2FoldChange > 0) |>
    select(stim, interval) |>
    {function(x) split(x, x$stim)}() |>
    map(~pull(., interval))


# Plots
n_peaks_plot <- 
    ggplot(peak_summary, aes(x = stim, y = value, fill = name)) +
    geom_col(position = "identity") +
    scale_fill_manual(values = c("#8856A7", "#E0E0E0")) +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 8)) +
    labs(x = NULL, y = "Number of peaks", fill = "N Peaks:",
         title = "Number of total peaks and tested peaks\n ")


library(eulerr)

fit <- euler(peak_lists)

euler_plot <- 
    plot(fit, quantities = list(cex = .6), labels = list(cex = .7), 
         fills = stim_colors$col, alpha = .5, edges = "grey50") |>
    wrap_elements(full = _, clip = FALSE) +
    ggtitle("Peak overlaps across conditions\n(vs Unstim 24hrs)") +
    theme(plot.title = element_text(size = 8))

plot_out <- 
    n_peaks_plot + plot_spacer() + euler_plot +
    plot_layout(widths = c(.6, .1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./supp_fig_atac_peakshare.png", plot_out, width = 6.5, height = 3.5)
