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

risk_var <- "rs7582694"

ld_vcf <- 
    "./data/chr2:190970120-192970120.vcf.gz" |>
    data.table::fread(skip = "#CHROM") |>
    as_tibble() |>
    unite("ID", c(ID, REF, ALT), sep = "-")

ld_plink <- 
    "./data/chr2:190970120-192970120_r2.ld" |>
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

langefeld <- 
    "./data/langefeld_stat4.tsv" |>
    data.table::fread() |>
    mutate(chrom = "2") |>
    select(chrom, pos, rsid, other_allele = alleleA, effect_allele = alleleB, p)

langefeld_grch38 <- 
    "./data/langefeld_hg38.bed" |>
    read_tsv(col_names = FALSE) |>
    select(chrom = X1, pos = X3, rsid = X4) |>
    mutate(chrom = str_remove(chrom, "chr")) |>
    distinct()

langefeld_ld <-
    langefeld |>
    left_join(ld_risk_var, join_by(rsid, other_allele == ref, effect_allele == alt)) |>
    group_by(chrom, pos, rsid) |>
    nest() |>
    ungroup() |>
    inner_join(langefeld_grch38, join_by(chrom, rsid)) |>
    select(chrom, pos = pos.y, rsid, data) |>
    unnest(cols = data) |>
    data.table::as.data.table()

loc_langefeld <- 
    locus(data = langefeld_ld, 
	  gene = 'STAT4',
	  flank = 0.75e5,
	  LD = "r2",
	  ens_db = ens_data)

loc_langefeld_ggplot <- 
    gg_scatter(loc_langefeld, 
	       labels = c("index", "rs11889341"),
	       nudge_y = .5, 
	       legend_pos = "topright") +
    labs(title = "Langefeld et al.")

# Gene tracks
g <- gg_genetracks(loc_langefeld) 

# plot
#dir.create("plots")

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
    GRanges("chr2", 
	    IRanges(min(loc_langefeld$data$pos), 
		    max(loc_langefeld$data$pos)))

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
    dplyr::filter(Chr == 2, padj <= 0.01) |>
    separate(contrast, c("stim1", "stim2"), sep = "vs") |>
    dplyr::filter(stim1 == "IL4_24" | stim2 == "IL4_24") |>
    mutate(stim = case_when(stim1 == "IL4_24" ~ stim2,
			    stim2 == "IL4_24" ~ stim1)) |>
    mutate(middle = (Start + End) / 2,
	   middle = round(middle)) |>
    select(stim, middle) |>
    inner_join(gr_df, join_by(stim, middle == bp)) |>
    mutate(lab = "*",
	   stim = factor(stim, levels = levels(gr_df$stim)))


atac <- 
    ggplot(gr_df) +
    geom_line(aes(x = bp, y = score, group = 1, color = stim),
	      linewidth = .5) +
    geom_text(data = da_peaks, 
	      aes(x = middle, y = score, label = "*"),
	      nudge_y = 0.1, size = 8, size.unit = "pt") +
    geom_vline(xintercept = 191079016, linetype = 2, linewidth = .25) +
    scale_x_continuous(limits = range(loc_langefeld$data$pos), 
		       labels = function(x) round(x/1e6L, 2),
		       breaks = scales::pretty_breaks(6)) +
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
    coord_cartesian(ylim = c(0, .25)) +
    labs(x = NULL)

ggsave("./plots/locuszoom.png", 
       loc_langefeld_ggplot / atac / g + plot_layout(heights = c(1, 1, 1)),
       width = 6, height = 7)

# Susie
susie <- 
    read_tsv("data/susie_stat4.tsv") |>
    left_join(langefeld_ld, join_by(rsid, ref == other_allele, alt == effect_allele))

gwas_stat4_plot <- 
    ggplot() +
    geom_point(data = langefeld_ld, 
	       aes(pos, -log10(p))) +
    geom_point(data = filter(susie, !is.na(cs)), 
	       aes(x = pos, -log10(p), color = cs),
	       shape = 21, size = 4, fill = NA, stroke = 1) +
    scale_x_continuous(labels = function(x) x/1e6L) +
    theme_minimal() +
    theme(legend.position = "top",
	  panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    coord_cartesian(xlim = range(loc_langefeld$data$pos)) +
    labs(x = "Chromosome 2 (Mb)",
	 y = expression("-log"["10"]("p")),
	 color = "Credible set:")

pip_plot <-
    ggplot(susie, aes(x = pos, y = pip)) +
    geom_point() +
    geom_point(data = filter(susie, !is.na(cs)), 
	       aes(x = pos, y = pip, color = cs),
	       shape = 21, size = 4, fill = NA, stroke = 1) +
    scale_x_continuous(labels = function(x) x/1e6L) +
    theme_minimal() +
    theme(legend.position = "none",
	  panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    coord_cartesian(xlim = range(loc_langefeld$data$pos)) +
    labs(x = "Chromosome 2 (Mb)",
	 y = "PIP")


ggsave("./plots/gwas_stat4.png", 
       gwas_stat4_plot / pip_plot / g, 
       width = 6, height = 8)





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

