library(tidyverse)
library(trackplot)

stim_colors <- 
    read_tsv("../figure_colors.txt", col_names = c("stim", "color")) |>
    mutate(stim = str_replace(stim, "_", " "),
	   stim = paste0(stim, "h")) |>
    deframe()


bigwigs <- 
    "./results/bwa/merged_replicate/bigwig" |>
    list.files(full.name = TRUE, pattern = "*.bigWig")

names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", basename(bigwigs))

bigwigs <- bigwigs[c("unst_0", "unst_24", "IL4_24", "TLR7_24", "BCR_24", "DN2_24")]

bigwigs <- 
    read_coldata(bws = bigwigs, build = "hg38") |>
    mutate(bw_sample_names = str_replace(bw_sample_names, "_", " "),
	   bw_sample_names = str_replace(bw_sample_names, "unst", "Unstim"),
	   bw_sample_names = paste0(bw_sample_names, "h"))

il12a_chrom <- "chr3"
il12a_start <- 159988835 - 1e4
il12a_end <- 159996019 + 5e4
il12a_locus <- glue::glue("{il12a_chrom}:{il12a_start}-{il12a_end}")

il12a_track <- track_extract(colData = bigwigs, loci = il12a_locus)

pdf("./il12a.pdf", width = 7, height = 5)
track_plot(summary_list = il12a_track, 
	   groupAutoScale = TRUE,
	   regions = data.frame(chr = "chr3", start = 160028076 - 1, end = 160028076),
	   boxcolalpha = 1, boxcol = "black",
	   cytoband_track_height = .5,
	   scale_track_height = 2,
	   col = stim_colors[bigwigs$bw_sample_names])
dev.off()

ikk_chrom <- "chr1"
ikk_pos <- 206470429
ikk_start <- ikk_pos - 1e4 
ikk_end <- ikk_pos + 1e4
ikk_locus <- glue::glue("{ikk_chrom}:{ikk_start}-{ikk_end}")

ikk_track <- track_extract(colData = bigwigs, loci = ikk_locus)

pdf("./ikkbe.pdf", width = 7, height = 5)
track_plot(summary_list = ikk_track, 
	   groupAutoScale = TRUE,
	   regions = data.frame(chr = ikk_chrom, start = ikk_pos - 1, end = ikk_pos),
	   boxcolalpha = 1, boxcol = "black",
	   cytoband_track_height = .5,
	   scale_track_height = 2,
	   col = stim_colors[bigwigs$bw_sample_names])
dev.off()


