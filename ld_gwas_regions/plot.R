library(tidyverse)
library(extrafont)
library(glue)
library(AnnotationHub)
library(locuszoomr)
library(patchwork)
library(rtracklayer)
library(ggrepel)

select <- dplyr::select
filter <- dplyr::filter
slice <- dplyr::slice

# Colors
stim_colors <- 
    c("unst_0" = "grey80", 
      "unst_24" = "grey50",
      "IL4_24" = "black", 
      "BCR_24" = "#00a0e4",
      "TLR7_24" = "#488f31",
      "DN2_24" = "#de425b"
      )


ah <- AnnotationHub()
#query(ah, "EnsDb.Hsapiens.v105")

ens_data <- ah[["AH98047"]]

langefeld_regions <- 
    "./data/langefeld_regions.tsv" |>
    read_tsv(col_names = c("region", "gene"))

langefeld <- 
    "./data/langefeld_summ_stats.tsv" |>
    data.table::fread() |>
    select(gene_region, chr, pos, rsid, other_allele = alleleA, effect_allele = alleleB, p)

langefeld_hg38 <- 
    "./data/langefeld_hg38.bed" |>
    read_tsv(col_names = FALSE) |>
    select(chr = X1, pos = X3, rsid = X4) |>
    mutate(chr = str_remove(chr, "chr")) |>
    distinct()

sentinels <- 
    "./data/langefeld_sentinels.tsv" |>
    read_tsv()

# Region 17
langefeld_stats <- 
    langefeld |>
    filter(grepl("DEF6", gene_region)) |>
    select(-gene_region)

risk_var <- 
    sentinels |> 
    filter(grepl("IL12A", gene_region)) |>
    pull(snp_id)

region <-
    langefeld_regions |>
    filter(grepl("IL12A", gene)) |>
    pull(region)

ld_vcf <- 
    glue("./data/chr{region}.vcf.gz") |>
    data.table::fread(skip = "#CHROM") |>
    as_tibble() |>
    unite("ID", c(ID, REF, ALT), sep = "-")

ld_plink <- 
    glue("./data/chr{region}_r2.ld") |>
    data.table::fread() |>
    as_tibble() |>
    setNames(ld_vcf$ID) |>
    add_column(snp_id = ld_vcf$ID, .before = 1)

ld_risk_var <-
    ld_plink |>
    filter(grepl(risk_var, snp_id)) |>
    pivot_longer(-snp_id, names_to = "var_id", values_to = "r2") |>
    select(var_id, r2) |>
    separate(var_id, c("rsid", "ref", "alt"), sep = "-", convert = TRUE)

langefeld_ld <-
    langefeld_stats |>
    mutate(chr = as.character(chr)) |>
    left_join(ld_risk_var, join_by(rsid, other_allele == ref, effect_allele == alt)) |>
    group_by(chr, pos, rsid) |>
    nest() |>
    ungroup() |>
    inner_join(langefeld_hg38, join_by(chr, rsid)) |>
    select(chr, pos = pos.y, rsid, data) |>
    unnest(cols = data) |>
    data.table::as.data.table()
    
locus_i <- genes(ens_data, filter = GeneNameFilter("IL12A"))

flank_left <- start(locus_i) - min(langefeld_ld$pos)
flank_right <- max(langefeld_ld$pos) - end(locus_i)

loc <- 
    locus(data = langefeld_ld, 
	  gene = 'IL12A',
	  flank = c(flank_left + 1, flank_right + 1) - 2e5,
	  LD = "r2",
	  ens_db = ens_data)


pip_df <- 
    "data/susie_IL12A_pip.tsv" |>
    read_tsv(col_types = "icccdcc")

plot_df <-
    left_join(langefeld_ld, pip_df,
	      join_by(rsid, other_allele == ref, effect_allele == alt)) |>
    mutate(r2interval = case_when(r2 >= 0 & r2 < .2 ~ "0.0 - 0.2",
				  r2 >= .2 & r2 < .4 ~ "0.2 - 0.4",
				  r2 >= .4 & r2 < .6 ~ "0.4 - 0.6",
				  r2 >= .6 & r2 < .8 ~ "0.6 - 0.8",
				  r2 >= .8 & r2 <= 1 ~ "0.8 - 1.0")) |>
    arrange(r2) |>
    mutate(r2interval = fct_rev(r2interval))

gwas_plot <- 
    ggplot() +
    geom_point(data = plot_df, 
	       aes(x = pos, y = -log10(p), fill = r2interval), 
	       size = 2, shape = 21, stroke = .25) +
    geom_point(data = filter(plot_df, !is.na(cs)), 
	       aes(x = pos, -log10(p), color = cs),
	       shape = 21, size = 4, fill = NA, stroke = 1) +
    scale_x_continuous(limits = loc$xrange,
		       labels = function(x) x/1e6L) +
    scale_color_manual(values = c("black", "Blue Violet", "blue")) +
    scale_fill_manual(values = c("0.0 - 0.2" = "#486CD9",
				 "0.2 - 0.4" = "#6BEBEC",
				 "0.4 - 0.6" = "#5DC83B",
				 "0.6 - 0.8" = "#F3A83B",
				 "0.8 - 1.0" = "#EB3223")) +
    theme_minimal() +
    theme(axis.title = element_text(size = 10),
	  panel.grid = element_blank(),
	  legend.title = element_text(size = 9),
	  legend.text = element_text(size = 9),
	  legend.key.height = unit(.5, "lines"),
	  legend.spacing.y = unit(-0.25, "cm"),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = NULL,
	 y = expression("-log"["10"]("p")),
	 color = "Credible set:",
	 fill = "r2:")

pip_plot <-
    ggplot() +
    geom_point(data = plot_df,
	       aes(x = pos, y = pip, fill = r2interval),
	       size = 2, shape = 21, stroke = .25) +
    geom_point(data = filter(plot_df, !is.na(cs)), 
	       aes(x = pos, y = pip, color = cs),
	       shape = 21, size = 4, fill = NA, stroke = 1) +
    geom_label_repel(data = plot_df |> filter(!is.na(cs)) |> group_by(cs) |> top_n(1, pip) |> ungroup(),
		    aes(x = pos, y = pip, label = rsid),
		    min.segment.length = 0, segment.size = .25, 
		    size = 3, alpha = .75) +
    scale_x_continuous(limits = loc$xrange,
		       labels = function(x) x/1e6L) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    scale_color_manual(values = c("black", "Blue Violet", "blue")) +
    scale_fill_manual(values = c("0.0 - 0.2" = "#486CD9",
				 "0.2 - 0.4" = "#6BEBEC",
				 "0.4 - 0.6" = "#5DC83B",
				 "0.6 - 0.8" = "#F3A83B",
				 "0.8 - 1.0" = "#EB3223")) +
    theme_minimal() +
    theme(axis.title = element_text(size = 10),
	  legend.position = "none",
	  panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = NULL,
	 y = "PIP") +
    coord_cartesian(clip = "off")

# ATAC-seq
bigwigs <- 
    "../atacseq/results/bwa/merged_replicate/bigwig" |>
    list.files(full.name = TRUE, pattern = "*.bigWig")

names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", basename(bigwigs))

da_files <- 
    list.files("../atacseq/results_deseq2", 
	       pattern = "\\.deseq2\\.FDR0\\.05\\.results\\.txt",
	       full.names = TRUE)

da_names <- 
    basename(da_files) |>
    strsplit("\\.") |>
    map_chr(1)

da_all <- da_files |>
    setNames(da_names) |>
    {function(x) x[!grepl("72", x)]}() |>
    map_dfr(~read_tsv(.) |> select(Geneid:padj), .id = "contrast") |>
    dplyr::filter(grepl("^chr", Chr)) |>
    mutate(Chr = str_remove(Chr, "chr"))

interv <- 
    GRanges(paste0("chr", unique(ld_vcf[[1]])), 
	    IRanges(loc$xrange[1], loc$xrange[2])) 

gr <- 
    map(bigwigs, ~import(., which = interv)) |>
    map_dfr(~as.data.frame(.) |> as_tibble(), .id = "stim") |>
    mutate(stim = factor(stim, levels = names(stim_colors)))

gr_df <- 
    gr |>
    mutate(score = score/max(score),
	   bp = map2(start, end, ~.x:.y)) |>
    select(stim, bp, score) |>
    unnest(cols = bp)

da_peaks <-
    da_all |>
    filter(Chr == 3, padj <= 0.01) |>
    separate(contrast, c("stim1", "stim2"), sep = "vs") |>
    filter(stim1 == "IL4_24" | stim2 == "IL4_24") |>
    mutate(stim = case_when(stim1 == "IL4_24" ~ stim2,
			    stim2 == "IL4_24" ~ stim1)) |>
    mutate(middle = (Start + End) / 2,
	   middle = round(middle)) |>
    select(stim, middle) |>
    inner_join(gr_df, join_by(stim, middle == bp)) |>
    mutate(lab = "*",
	   stim = factor(stim, levels = levels(gr_df$stim)))

atac_plot <- 
    ggplot(gr_df) +
    geom_line(aes(x = bp, y = score, group = 1, color = stim),
	      linewidth = .75) +
    geom_text(data = da_peaks, 
	      aes(x = middle, y = score, label = "*"),
	      nudge_y = 0.1, size = 10, size.unit = "pt") +
    geom_vline(xintercept = langefeld_ld |> filter(rsid == risk_var) |> pull(pos), 
	       linetype = 2, linewidth = .25) +
    geom_vline(xintercept = langefeld_ld |> filter(rsid == risk_var) |> pull(pos), 
	       linetype = 2, linewidth = .25) +
    scale_x_continuous(limits = loc$xrange, 
		       labels = function(x) round(x/1e6L, 2)) +
    scale_color_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 1) +
    theme_minimal() +
    theme(axis.ticks.y = element_blank(),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  strip.text = element_blank(),
	  legend.position = "none",
	  panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    coord_cartesian(ylim = c(0, .5)) +
    labs(x = NULL)

g <- gg_genetracks(loc, cex.text = 0.7, cex.axis = 0.9)

p1 <- plot_spacer()
p2 <- (gwas_plot / pip_plot / atac_plot / g + plot_layout(heights = c(1, .5, 1, .5))) 
p3 <- plot_spacer()

p_out <- 
    p1 / p2 / p3 + 
    plot_annotation(title = "IL12A") +
    plot_layout(heights = c(.1, 1, .1))

ggsave("./plots/il12a.png", p_out, width = 8.5, height = 8.5)





# hg19
library(EnsDb.Hsapiens.v75)

sentinels <- 
    "./data/langefeld_sentinels.tsv" |>
    read_tsv()

regions_df <-
    "./data/langefeld_regions.tsv" |>
    read_tsv(col_names = c("region", "gene_region"))

plot_susie <- function(ix) {
    
    gene_ix <- regions_df$gene_region[ix]
    region_ix <- regions_df$region[ix]

    risk_var <- 
	sentinels |>
	filter(gene_region == gene_ix) |>
	pull(snp_id)

    ld_vcf <- 
	glue("./data/chr{region_ix}.vcf.gz") |>
	data.table::fread(skip = "#CHROM") |>
	as_tibble() |>
	unite("ID", c(ID, REF, ALT), sep = "-")

    ld_plink <- 
	glue("./data/chr{region_ix}_r2.ld") |>
	data.table::fread() |>
	as_tibble() |>
	setNames(ld_vcf$ID) |>
	add_column(snp_id = ld_vcf$ID, .before = 1)

    ld_risk_var <-
	ld_plink |>
	filter(grepl(paste0(risk_var, "-"), snp_id)) |>
	pivot_longer(-snp_id, names_to = "var_id", values_to = "r2") |>
	select(var_id, r2) |>
	separate(var_id, c("rsid", "ref", "alt"), sep = "-", convert = TRUE)

    stats <- 
	read_tsv("./data/langefeld_summ_stats.tsv") |>
	filter(gene_region == gene_ix) |>
	select(chr, pos, rsid, alleleA, alleleB, p) |>
	left_join(ld_risk_var, join_by(rsid, alleleA == ref, alleleB == alt)) |>
	data.table::as.data.table()

    gene_lab <- str_split(gene_ix, "-")[[1]][1] 
    
    edb <- get("EnsDb.Hsapiens.v75")
    locus_i <- genes(edb, filter = GeneNameFilter(gene_lab))
    
    flank_left <- start(locus_i) - min(stats$pos)
    flank_right <- max(stats$pos) - end(locus_i)

    loc <- 
	locus(data = stats, 
	      gene = gene_lab,
	      flank = c(flank_left + 1, flank_right + 1),
	      LD = "r2",
	      ens_db = "EnsDb.Hsapiens.v75")

    # Gene tracks
    g <- gg_genetracks(loc, cex.text = 0.7, cex.axis = 0.9)

    pip_df <- 
	glue("data/susie_{gene_ix}_pip.tsv") |>
	read_tsv(col_types = "icccdcc") |>
	left_join(stats, join_by(rsid, ref == alleleA, alt == alleleB)) |>
	mutate(r2interval = case_when(r2 >= 0 & r2 < .2 ~ "0.0 - 0.2",
				      r2 >= .2 & r2 < .4 ~ "0.2 - 0.4",
				      r2 >= .4 & r2 < .6 ~ "0.4 - 0.6",
				      r2 >= .6 & r2 < .8 ~ "0.6 - 0.8",
				      r2 >= .8 & r2 <= 1 ~ "0.8 - 1.0")) |>
	arrange(r2) |>
	mutate(r2interval = fct_rev(r2interval))

    gwas_plot <- 
	ggplot() +
	geom_point(data = pip_df, 
		   aes(x = pos, y = -log10(p), fill = r2interval), 
		   size = 2, shape = 21, stroke = .25) +
	geom_point(data = filter(pip_df, !is.na(cs)), 
		   aes(x = pos, -log10(p), color = cs),
		   shape = 21, size = 4, fill = NA, stroke = 1) +
	scale_x_continuous(limits = loc$xrange,
			   labels = function(x) x/1e6L) +
	scale_color_manual(values = c("black", "Blue Violet", "blue")) +
	scale_fill_manual(values = c("0.0 - 0.2" = "#486CD9",
				     "0.2 - 0.4" = "#6BEBEC",
				     "0.4 - 0.6" = "#5DC83B",
				     "0.6 - 0.8" = "#F3A83B",
				     "0.8 - 1.0" = "#EB3223")) +
	theme_minimal() +
	theme(axis.title = element_text(size = 10),
	      panel.grid = element_blank(),
	      legend.title = element_text(size = 9),
	      legend.text = element_text(size = 9),
	      legend.key.height = unit(.5, "lines"),
	      legend.spacing.y = unit(-0.25, "cm"),
	      plot.background = element_rect(fill = "white", color = "white")) +
	labs(x = NULL,
	     y = expression("-log"["10"]("p")),
	     color = "Credible set:",
	     fill = "r2:")

    pip_plot <-
	ggplot() +
	geom_point(data = pip_df,
		   aes(x = pos, y = pip, fill = r2interval),
		   size = 2, shape = 21, stroke = .25) +
	geom_point(data = filter(pip_df, !is.na(cs)), 
		   aes(x = pos, y = pip, color = cs),
		   shape = 21, size = 4, fill = NA, stroke = 1) +
	geom_label_repel(data = pip_df |> filter(!is.na(cs)) |> group_by(cs) |> top_n(1, pip) |> ungroup(),
			aes(x = pos, y = pip, label = rsid),
			min.segment.length = 0, segment.size = .25, 
			size = 3, alpha = .75) +
	scale_x_continuous(limits = loc$xrange,
			   labels = function(x) x/1e6L) +
	scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
	scale_color_manual(values = c("black", "Blue Violet", "blue")) +
	scale_fill_manual(values = c("0.0 - 0.2" = "#486CD9",
				     "0.2 - 0.4" = "#6BEBEC",
				     "0.4 - 0.6" = "#5DC83B",
				     "0.6 - 0.8" = "#F3A83B",
				     "0.8 - 1.0" = "#EB3223")) +
	theme_minimal() +
	theme(axis.title = element_text(size = 10),
	      legend.position = "none",
	      panel.grid = element_blank(),
	      plot.background = element_rect(fill = "white", color = "white")) +
	labs(x = NULL,
	     y = "PIP")
    
    condz_df <- 
	glue("./data/susie_{gene_ix}_krss.tsv") |>
	read_tsv() |>
	mutate(flag = logLR > 2 & abs(z) > 2)

    condz_plot <- 
        ggplot(condz_df, aes(x = condmean, y = z, color = flag)) +
        geom_abline() +
        geom_point(size = .5) +
	scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
        theme_bw() +
	theme(legend.position = "none") +
        labs(x = "Expected value", y = "Observed value",
	     title = "Susie 'kriging_rss'") +
	coord_fixed()

    p1 <- plot_spacer()
    p2 <- (gwas_plot / pip_plot / g) 
    p3 <- plot_spacer()
    p4 <- (condz_plot + plot_spacer() + plot_spacer())
    p5 <- plot_spacer()
    
    p_out <- 
	p1 / p2 / p3 / p4 + 
	plot_annotation(title = gene_ix) +
	plot_layout(heights = c(.1, 1, .1, .5, .1))
    
    ggsave(glue("./plots/susie_{ix}.pdf"), p_out, width = 8.5, height = 11)
}

walk(1:nrow(regions_df), plot_susie)


# IL12A #####################################################################
stim_colors <- 
    read_tsv("../figure_colors.txt", col_names = c("stim", "color")) |>
    mutate(stim = str_replace(stim, "_", " "),
	   stim = paste0(stim, "h")) |>
    deframe()

summ_stats <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/harmonized/26502338-GCST003156-EFO_0002690.h.tsv.gz" |>
    data.table::fread() |>
    select(chr = chromosome, pos = base_pair_location, rsid = variant_id, p = p_value)

risk_var <- 
    summ_stats |>
    filter(rsid == "rs485499") |>
    pull(pos)

loc <- 
    locus(data = summ_stats, 
	  gene = 'IL12A',
	  flank = c(1e5, 1.25e5),
	  ens_db = ens_data)

bigwigs <- 
    "../atacseq/results/bwa/merged_replicate/bigwig" |>
    list.files(full.name = TRUE, pattern = "*.bigWig")

names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", basename(bigwigs))

da_files <- 
    list.files("../atacseq/results_deseq2", 
	       pattern = "*vsunst_24\\.tsv$",
	       full.names = TRUE)

da_names <- 
    basename(da_files) |>
    strsplit("vs") |>
    map_chr(1)

da_all <- 
    da_files |>
    setNames(da_names) |>
    map_dfr(read_tsv, .id = "stim") |>
    filter(grepl("^chr", Chr)) |>
    mutate(Chr = str_remove(Chr, "chr"))

interv <- 
    GRanges("chr3", IRanges(loc$xrange[1], loc$xrange[2])) 

gr <- 
    map(bigwigs, ~import(., which = interv)) |>
    map_dfr(~as.data.frame(.) |> as_tibble(), .id = "stim") |>
    mutate(stim = str_replace(stim, "_", " "),
	   stim = str_replace(stim, "unst", "Unstim"),
	   stim = paste0(stim, "h"),
	   stim = factor(stim, levels = names(stim_colors)))

gr_df <- 
    gr |>
    mutate(score = score/max(score),
	   bp = map2(start, end, ~.x:.y)) |>
    select(stim, bp, score) |>
    unnest(cols = bp)

da_peaks <-
    da_all |>
    filter(Chr == 3, padj <= 0.01) |>
    mutate(middle = (Start + End) / 2,
	   middle = round(middle)) |>
    select(stim, middle) |>
    mutate(stim = str_replace(stim, "_", " "),
	   stim = str_replace(stim, "unst", "Unstim"),
	   stim = paste0(stim, "h"),
	   stim = factor(stim, levels = names(stim_colors))) |>
    inner_join(gr_df, join_by(stim, middle == bp))

atac_plot <- 
    ggplot(gr_df) +
    geom_line(aes(x = bp, y = score, group = 1, color = stim),
	      linewidth = .75) +
    geom_text(data = da_peaks, 
	      aes(x = middle, y = score, label = "*"),
	      nudge_y = 0.1, size = 10, size.unit = "pt") +
    geom_vline(xintercept = risk_var, linetype = 2, linewidth = .25) +
    scale_x_continuous(limits = loc$xrange, 
		       labels = function(x) round(x/1e6L, 2)) +
    scale_color_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 1) +
    theme_minimal() +
    theme(axis.ticks.y = element_blank(),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  strip.text = element_blank(),
	  legend.position = "none",
	  panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    coord_cartesian(ylim = c(0, .5)) +
    labs(x = NULL)

g <- gg_genetracks(loc, cex.text = 0.7, cex.axis = 0.9)

p_out <- 
    atac_plot / plot_spacer() / g + 
    plot_annotation(title = "IL12A") +
    plot_layout(heights = c(1, .1, .25))

ggsave("./plots/IL12A.png", p_out, height = 4)











# IL12A v2 #####################################################################
stim_colors <- 
    read_tsv("../figure_colors.txt", col_names = c("stim", "color")) |>
    mutate(stim = str_replace(stim, "_", " "),
	   stim = paste0(stim, "h")) |>
    deframe()

summ_stats <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/harmonized/26502338-GCST003156-EFO_0002690.h.tsv.gz" |>
    data.table::fread() |>
    select(chr = chromosome, pos = base_pair_location, rsid = variant_id, p = p_value)

loc <- 
    locus(data = summ_stats, 
	  gene = 'IL12A',
	  flank = c(1e4, .5e5),
	  ens_db = ens_data)

bigwigs <- 
    "../atacseq/results/bwa/merged_replicate/bigwig" |>
    list.files(full.name = TRUE, pattern = "*.bigWig")

names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", basename(bigwigs))

da_files <- 
    list.files("../atacseq/results_deseq2", 
	       pattern = "*vsunst_24\\.tsv$",
	       full.names = TRUE)

da_names <- 
    basename(da_files) |>
    strsplit("vs") |>
    map_chr(1)

da_all <- 
    da_files |>
    setNames(da_names) |>
    map_dfr(read_tsv, .id = "stim") |>
    filter(grepl("^chr", Chr)) |>
    mutate(Chr = str_remove(Chr, "chr"))

interv <- 
    GRanges(glue("chr{loc$seqname}"), IRanges(loc$xrange[1], loc$xrange[2])) 

gr <- 
    map(bigwigs, ~import(., which = interv)) |>
    map_dfr(~as.data.frame(.) |> as_tibble(), .id = "stim") |>
    mutate(stim = str_replace(stim, "_", " "),
	   stim = str_replace(stim, "unst", "Unstim"),
	   stim = paste0(stim, "h"),
	   stim = factor(stim, levels = names(stim_colors)))

gr_df <- 
    gr |>
    mutate(score = score/max(score),
	   bp = map2(start, end, ~.x:.y)) |>
    select(stim, bp, score) |>
    unnest(cols = bp)

da_peaks <-
    da_all |>
    filter(Chr == loc$seqname, padj <= 0.01) |>
    select(stim, interval, Chr, Start, End) |>
    mutate(stim = str_replace(stim, "_", " "),
	   stim = str_replace(stim, "unst", "Unstim"),
	   stim = paste0(stim, "h"),
	   stim = factor(stim, levels = names(stim_colors))) |>
    {function(x) inner_join(x = gr_df, y = x, join_by(stim, between(bp, Start, End)))}() |>
    group_by(stim, interval) |>
    slice_max(score) |>
    ungroup()

vars_df <-
    tribble(~rsid, ~pos,
	    "rs564799", 160011200,
	    "rs564976", 160011272, 
	    "rs485499", 160028076,
	    "rs55843920", 160011393
	    )

atac_plot <- 
    ggplot(gr_df) +
    geom_vline(xintercept = vars_df$pos, linetype = 2, linewidth = .25) +
    geom_ribbon(aes(x = bp, ymin = 0, ymax = score, fill = stim), 
		color = 'black', linewidth = .1) +
    geom_text(data = da_peaks, 
	      aes(x = bp, y = score, label = "*"),
	      size = 8, fontface = "bold", size.unit = "pt") +
    scale_x_continuous(limits = loc$xrange, 
		       labels = function(x) round(x/1e6L, 2)) +
    scale_y_continuous(expand = expansion(mult = c(0, .2))) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 1) +
    theme_minimal() +
    theme(axis.ticks.y = element_blank(),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  strip.text = element_blank(),
	  legend.position = "none",
	  panel.grid = element_blank(),
	  panel.spacing.y = unit(0, "lines"),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = NULL)

g <- gg_genetracks(loc, cex.text = 0.7, cex.axis = 0.9)

p_out <- 
    atac_plot / g + 
    plot_annotation(title = "IL12A") +
    plot_layout(heights = c(1, .25))

ggsave("./plots/IL12A.png", p_out, height = 3.5)



