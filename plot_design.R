library(tidyverse)
library(readxl)
library(patchwork)
library(ggridges)

stim_colors <- 
    c("unstim_0h" = "#dbdbdb", 
      "unstim_4h" = "#979797", 
      "unstim_24h" = "#777777",
      "IL4_4h" = "#585858",
      "IL4_24h" = "#3b3b3b",
      "IL4_72h"  = "#000000",
      "sCD40L_4h" = "#ffeba6",
      "sCD40L_24h" = "#ffde68",
      "sCD40L_48h" = "goldenrod3",
      "sCD40L_72h" = "goldenrod4",
      "BCR_4h" = "#a1c2ed",
      "BCR_24h" = "#6996e3",
      "BCR_48h" = "#4060c8",
      "BCR_72h" = "#1a318b",
      "TLR7_4h" = "#98ab76",
      "TLR7_24h" = "#748f46",
      "TLR7_48h" = "#47632a",
      "TLR7_72h" = "#275024",
      "TLR9_4h" = "#b695bc",
      "TLR9_24h" = "#90719f",
      "TLR9_48h" = "#574571",
      "TLR9_72h" = "#2d223c",
      "BCR+TLR7_4h" = "#ffa7db",
      "BCR+TLR7_24h" = "#ff86d0",
      "BCR+TLR7_48h" = "#ff61c4",
      "BCR+TLR7_72h" = "#ff2bb8",
      "DN2_4h" = "#e6907a",
      "DN2_24h" = "#d76b51",
      "DN2_48h" = "#c5432a",
      "DN2_72h" = "#b00000")

lowinput <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_lowinput/metadata/description_hsapiens.tsv" |>
    read_tsv() |>
    distinct(donor_id = sample_id, stim, time) |>
    mutate(donor_id = sub("MGB-", "", donor_id),
	   donor_id = sub("\\.rep\\.\\d", "", donor_id),
	   time = sub("rs", "", time)) |>
    arrange(donor_id, stim, time) |>
    distinct() |>
    filter(donor_id != "BLANK")

rnaseq <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/metadata.tsv" |>
    read_tsv(col_types = c(.default = "c")) |>
    distinct(donor_id = subject_id, stim) |>
    mutate(stim = recode(stim, "unstday0" = "unstim"),
	   time = ifelse(stim == "unstim", "0h", "24h"))

citeseq <- 
    "./data/Lupus.xlsx" |>
    read_excel(sheet = 4, skip = 1, col_names = FALSE, col_types = "text") |>
    fill(1) |>
    pivot_longer(-1) |>
    drop_na(value) |>
    select(donor_id = 1, value) |>
    separate(value, c("stim", "time"), sep = " ") |>
    distinct()

atac <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/atacseq/samplesheet.csv" |>
    read_csv() |>
    distinct(fastq_1) |>
    mutate(fastq_1 = basename(fastq_1)) |>
    extract(fastq_1, c("donor_id", "stim"), "\\d+_([^_]+)_([^_]+_[^_]+)_.*") |>
    separate(stim, c("stim", "time"), sep = "_") |>
    mutate(time = recode(time, "day0" = "0h"),
	   stim = recode(stim, "unst" = "unstim")) |>
    filter(donor_id != "3donors") |>
    arrange(donor_id, stim, time)

assay_df <- 
    bind_rows("Low-input RNAseq" = lowinput,
	      "Normal-input RNAseq" = rnaseq,
	      "CITEseq" = citeseq,
	      "ATACseq" = atac, 
	      .id = "assay") |>
    unite("color", c(stim, time), sep = "_", remove = FALSE) |>
    filter(!(time == 0 & stim != "unstim")) |>
    mutate(color = factor(color, levels = names(stim_colors)),
	   time = factor(parse_number(time), levels = c(0, 4, 24, 48, 72)),
	   stim = factor(stim, levels = c("unstim", "IL4", "sCD40L", "TLR9", "TLR7", "BCR", "BCR+TLR7", "DN2")),
	   assay = factor(assay, levels = c("Low-input RNAseq", "Normal-input RNAseq", "ATACseq", "CITEseq")))

summ_assays <- distinct(assay_df, assay, color, stim, time)

summ_assays <- 
    bind_rows(summ_assays,
	      distinct(summ_assays, assay) |>
	      mutate(color = NA, stim = factor("IL4"), time = factor(48)))

# Plot

sample_size_p <- 
    assay_df |>
    distinct(assay, donor_id) |>
    count(assay) |>
    ggplot(aes(x = n, y = assay)) +
    geom_col(fill = "black") +
    scale_x_reverse(limits = c(NA, 0)) +
    geom_text(aes(label = n), hjust = 0, color = "white", fontface = "bold") +
    facet_grid(assay~., switch = "y", scales = "free") +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  panel.background = element_rect(color = "white", fill = "white"),
	  axis.title = element_blank(),
	  axis.ticks = element_blank(),
	  axis.text = element_blank(),
	  strip.text.y.left = element_text(angle = 0, size = 12),
	  strip.placement = "outside",
	  strip.background = element_rect(fill = "white", color = "white"),
	  panel.spacing = unit(0, "lines"),
	  panel.border = element_blank(),
	  plot.margin = margin(.5, 0, .5, .5))

p0 <- 
    tibble(assay = fct_inorder(levels(assay_df$assay)), g = c(0, 1, 1, 2)) |>
    ggplot(aes(x = 1, y = assay, group = g)) +
	geom_point() +
	geom_line() +
	theme_minimal() +
	theme(panel.grid = element_blank(),
	      panel.background = element_rect(color = "white", fill = "white"),
	      axis.title = element_blank(),
	      axis.ticks = element_blank(),
	      axis.text = element_blank(),
	      plot.margin = margin(.5, 0, .5, 0))

p <- 
    ggplot(summ_assays, aes(x = time, y = assay, fill = color)) +
    geom_tile(color = NA) +
    scale_fill_manual(values = stim_colors, na.value = "white") +
    scale_x_discrete(position = "top") +
    facet_grid(assay~stim, switch = "y", scales = "free", space = "free") +
    theme_bw() +
    theme(legend.position = "none",
	  panel.grid = element_blank(),
	  panel.background = element_rect(color = "white", fill = "white"),
	  axis.title = element_blank(),
	  axis.text.y = element_blank(),
	  axis.ticks = element_blank(),
	  strip.text.x = element_text(size = 12),
	  strip.text.y.left = element_blank(),
	  strip.placement = "outside",
	  panel.spacing = unit(0, "lines"),
	  panel.border = element_rect(linewidth = .1),
	  strip.background = element_rect(fill = "white", color = "white"),
	  plot.margin = margin(.5, .5, .5, 0))

ggsave("./study_design.png", 
       sample_size_p + p0 + p + plot_layout(widths = c(.1, .025, 1)), 
       height = 1.8, width = 10, dpi = 600)

# density
#y <- rnorm(5e4, 100, 10)
#
#cutoff <- quantile(y, probs = 0.95)
#
#hist_y <- density(y, from = 50, to = 150) |>
#    broom::tidy() |>
#    mutate(area = x >= cutoff)
#
hist_df <- tibble(x = rnorm(5e4, 100, 10), y = factor("dummy"))

p2 <- 
    ggplot(data = hist_df) +
    stat_density_ridges(aes(x = x, y = y, fill = factor(after_stat(quantile))),
			geom = "density_ridges_gradient",
			calc_ecdf = TRUE,
			quantiles = 0.9) + 
    theme_minimal() +
    scale_y_discrete(expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(values = c("1" = "grey80", "2" = "magenta4")) +
    theme(plot.background = element_rect(color = "white", fill = "white"),
	  axis.title.x = element_text(size = 10),
	  axis.title.y = element_blank(),
	  axis.text = element_blank(),
	  panel.grid = element_blank(),
	  legend.position = "none") +
    coord_cartesian(clip = "off") +
    labs(x = "Number of heterozygous risk variants per individual")

ggsave("./sledist.png", p2, height = 2, width = 4, dpi = 600)
