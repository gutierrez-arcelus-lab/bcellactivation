library(tidyverse)
library(patchwork)

stim_colors <- 
    c("unstday0" = "grey80", 
      "BCR" = "blue",
      "TLR7" = "#488f31",
      "DN2" = "#de425b"
      )

# GWAS regions
regions <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/colocalization/finemap/data/regions_bentham_1mb.tsv" |>
    read_tsv(col_names = c("region", "locus")) |>
    extract(region, c("chr", "left", "right"), 
	    "(chr[^:]+):(\\d+)-(\\d+)", 
	    convert = TRUE)

# Gene tracks
tracks <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/gencode_v38_gene_structure.tsv" |>
    read_tsv() |>
    group_by(gene_id) |>
    mutate(tss = ifelse(strand == "+", min(start), max(end))) |>
    ungroup() |>
    inner_join(regions, join_by(chr, between(tss, left, right))) |>
    select(chr, gwas_locus = locus, gene_id, gene_name, transcript_id, 
	   strand, feature, i, start, end)

# ASE in genes on GWAS regions
ase <- read_tsv("./ase_gwasregions.tsv", col_types = "cccdciiddccc") |>
    group_by(variantID) |>
    filter(any(q < 0.05)) |>
    ungroup()

itgam <- filter(ase, locus == "ITGAM")


itgax_tracks <- filter(tracks, gene_name == "ITGAX")
itgam_itgax <- filter(itgam, gene_name == "ITGAX") |>
    mutate(s = q < 0.05)

p <- 
    ggplot(data = itgam_itgax,
	   aes(x = position, y = refCount/total)) +
    geom_jitter(aes(color = stim, fill = stim, shape = s), 
	       size = 4, width = 100) +
    scale_x_continuous(limits = c(min(itgax_tracks$start), max(itgax_tracks$end)),
		       breaks = unique(itgam_itgax$position),
		       labels = function(x) round(x/1e6L, 1)) +
    scale_y_continuous(limits = c(-0.1, 1.1), 
		       breaks = c(0, .5, 1)) +
    scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 1), guide = "none") +
    scale_color_manual(values = stim_colors) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(sample_id~., ncol = 1, strip.position = "right") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(), 
	  panel.grid.minor.y = element_blank(),
	  axis.text.x = element_text(angle = 90, hjust = 0, vjust = 1),
	  strip.text = element_text(size = 9),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = NULL)

t <- 
    ggplot(data = itgax_tracks) +
    geom_segment(data = filter(itgax_tracks, feature == "exon"),
		 aes(x = start, xend = end, y = gene_name, yend = gene_name, group = i),
		 linewidth = 10) +
    geom_segment(data = filter(itgax_tracks, feature == "intron"),
		 aes(x = start, xend = end, y = gene_name, yend = gene_name, group = i),
		 linewidth = 1) +
    scale_x_continuous(limits = c(min(itgax_tracks$start), max(itgax_tracks$end)),
		       breaks = unique(itgam_itgax$position)) +
    theme_minimal() +
    theme(panel.grid.minor.x = element_blank(), 
	  panel.grid.major.y = element_blank(),
	  axis.text.x = element_blank(),
	  axis.title = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white"))

out <- p / t + plot_layout(heights = c(1, .25))
ggsave("./plots/itgax.png", out, width = 14, height = 12)


itgam_tracks <- filter(tracks, gene_name == "ITGAM")
itgam_itgam <- filter(itgam, gene_name == "ITGAM") |>
    mutate(s = q < 0.05)

p <- 
    ggplot(data = itgam_itgam,
	   aes(x = position, y = refCount/total)) +
    geom_jitter(aes(color = stim, fill = stim, shape = s), 
	       size = 4, width = 100) +
    scale_x_continuous(limits = c(min(itgam_tracks$start), max(itgam_tracks$end)),
		       breaks = unique(itgam_itgam$position),
		       labels = function(x) round(x/1e6L, 1)) +
    scale_y_continuous(limits = c(-0.1, 1.1), 
		       breaks = c(0, .5, 1)) +
    scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 1), guide = "none") +
    scale_color_manual(values = stim_colors) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(sample_id~., ncol = 1, strip.position = "right") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(), 
	  panel.grid.minor.y = element_blank(),
	  axis.text.x = element_text(angle = 90, hjust = 0, vjust = 1),
	  strip.text = element_text(size = 9),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = NULL)

t <- 
    ggplot(data = itgam_tracks) +
    geom_segment(data = filter(itgam_tracks, feature == "exon"),
		 aes(x = start, xend = end, y = gene_name, yend = gene_name, group = i),
		 linewidth = 10) +
    geom_segment(data = filter(itgam_tracks, feature == "intron"),
		 aes(x = start, xend = end, y = gene_name, yend = gene_name, group = i),
		 linewidth = 1) +
    scale_x_continuous(limits = c(min(itgam_tracks$start), max(itgam_tracks$end)),
		       breaks = unique(itgam_itgam$position)) +
    theme_minimal() +
    theme(panel.grid.minor.x = element_blank(), 
	  panel.grid.major.y = element_blank(),
	  axis.text.x = element_blank(),
	  axis.title = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white"))

out <- p / t + plot_layout(heights = c(1, .25))
ggsave("./plots/itgam.png", out, width = 14, height = 12)

# Check if genotypes at risk variant associates with ASE in the region
fix_samples <- read_tsv("./mbv/fixed_ase/data/fix_metadata.tsv", col_names = FALSE) |>
    extract(X4, "wrong_id", "\\d+_(\\d+)_.+") |>
    select(donor_id = X1, wrong_id) |>
    distinct()

riskvar <- "chr16:31265490:G:A" 

vcf <- "./data/allchr.r2filtered.mgb.vcf.gz" |>
    read_tsv(comment = "##")

vcf2 <- "./mbv/fixed_ase/data/allchr.r2filtered.mgb.vcf.gz" |>
    read_tsv(comment = "##")

vcf2 |> 
    filter(`#CHROM` == "chr16", POS == 31265490)

vcf2_ori <- 
    list.files("/temp_work/ch229163/vcf/mbv/", 
	       pattern = "chr16\\.MGB\\.040\\d\\.vcf\\.gz$",
	       full.names = TRUE) |>
    map_dfr(~read_tsv(., comment = "##") |> 
	    filter(ID == riskvar) |>
	    pivot_longer(-(1:9), names_to = "donor_id")) |>
    extract(donor_id, c("donor_id"), "\\d+_[^-]+-(\\d+)") |>
    extract(FORMAT, "FORMAT", "([^:]+)") |>
    extract(value, "value", "([^:]+)")
    



vcf |> 
    filter(`#CHROM` == "chr16", POS == 31265490) |>
    pivot_longer(-(1:9)) |>
    extract(name, c("donor_id"), "\\d+_[^-]+-(\\d+)") |>
    filter(! donor_id %in% fix_samples$wrong_id) |>
    bind_rows(vcf2_ori)















