# Use 48Gb of RAM

Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)

library(tidyverse)
library(susieR)
library(cowplot)
library(gridExtra)
library(furrr)
library(ggrepel)

# Functions
run_susie <- function(locus_name) {

    region_id <- regions$region[regions$locus == locus_name] 

    vcf <- sprintf("./data/%s.adj.vcf.gz", region_id) |>
	read_tsv(comment = "##") 

    plink <- sprintf("./data/%s.ld", region_id) |>
	read_delim(delim = " ", col_names = FALSE)

    plink_filt <- plink[, sapply(plink, function(x) !all(is.na(x)))]
    ldmat <- data.matrix(plink_filt)
    rownames(ldmat) <- vcf$ID
    colnames(ldmat) <- vcf$ID

    summ_stats <- summ_stats_all |>
	filter(locus == locus_name) |>
	unite("ID", c(chr, bp, ref, alt), sep = ":") |>
	filter(ID %in% vcf$ID) |>
	mutate(z = beta/se)

    #condz <- kriging_rss(summ_stats$z, ldmat, n = sample_size, r_tol = 1e-04)
    #png(sprintf("./plots/diagnostic_%s.png", locus_name))
    #condz$plot
    #dev.off()

    fit <- susie_rss(summ_stats$z, 
		     ldmat, 
		     n = sample_size, 
		     L = 5, 
		     coverage = 0.9, 
		     r_tol = 1e-05)

    # if it does not converge, remove problematic SNPs and repeat until convergence
    iter <- 0L
    while ( !fit$converged & iter <= 20 ) {

	condz <- kriging_rss(summ_stats$z, ldmat, n = sample_size, r_tol = 1e-04)

	snp_rm <- as_tibble(condz$conditional_dist) |>
	    rowid_to_column() |>
	    arrange(desc(abs(z_std_diff))) |>
	    slice(1) |>
	    pull(rowid)

	ldmat <- ldmat[-snp_rm, -snp_rm]
	summ_stats <- slice(summ_stats, -snp_rm)

	fit <- susie_rss(summ_stats$z, 
			 ldmat, 
			 n = sample_size, 
			 L = 5, 
			 coverage = 0.9, 
			 r_tol = 1e-05)

	iter <- iter + 1L
    }

    return(fit)
}

make_pip_df <- function(fit_obj) {

    if (! is.null(fit_obj$sets$coverage) ) {

	coverage_df <- 
	    tibble(coverage = scales::percent(round(fit_obj$sets$coverage, 3)),
		   cs = paste0("L", fit_obj$sets$cs_index))

	cs_df <- enframe(fit_obj$sets$cs, "cs", "rowid") |>
	    unnest(cols = rowid) |>
	    left_join(coverage_df, by = "cs")

	pip_df <- enframe(fit_obj$pip, "ID", "pip") |>
	    separate(ID, c("chr", "pos", "ref", "alt"), sep = ":", convert = TRUE) |>
	    rowid_to_column() |>
	    left_join(cs_df) |>
	    mutate(cs = ifelse(!is.na(cs), paste0(cs, " (", coverage, ")"), NA))

    } else {
	
	pip_df <- 
	    enframe(fit_obj$pip, "ID", "pip") |>
	    separate(ID, c("chr", "pos", "ref", "alt"), sep = ":", convert = TRUE) |>
	    mutate(cs = NA)
    }

    return(pip_df)
}

plot_susie <- function(res_df) {


    tracks_df <- tracks |>
	filter(gwas_locus == unique(res_df$locus)) |>
	select(chr, locus = gene_name, feature, i, start, end) |>
	group_by(chr, locus) |>
	nest() |>
	ungroup()

    tracks_df$rowid <- rep(1:5, length.out = nrow(tracks_df))
    tracks_df <- unnest(tracks_df, cols = data)
    
    strand_df <- distinct(tracks, gene_name, strand)
   
    arrows_df <- 
	left_join(tracks_df, strand_df, by = join_by(locus == gene_name)) |>
	group_by(locus) |>
	mutate(pos = ifelse(strand == "+", max(end), min(start)) ) |>
	distinct(locus, strand, pos) |>
	mutate(x1 = ifelse(strand == "+", pos + 1e3, pos - 1e3))
    
    genes_df <- 
	left_join(tracks_df, strand_df, by = join_by(locus == gene_name)) |>
	group_by(locus) |>
	mutate(pos = ifelse(strand == "-", max(end), min(start)) ) |>
	ungroup() |>
	distinct(locus, pos) |>
	left_join(distinct(tracks_df, locus, rowid), join_by(locus))

    max_pip <- res_df |> 
	filter(!is.na(cs), stat == "Susie PIP") |>
	group_by(cs) |>
	slice_max(value) |>
	mutate(stat = paste0(stat, ";GWAS p-value")) |>
	separate_rows(stat, sep = ";")


    pip_plot <- 
	ggplot(data = res_df, 
	       aes(x = pos, y = value)) +
	#geom_vline(data = max_pip, aes(xintercept = pos), linetype = 2, color = "grey35") +
	geom_point() +
	geom_point(data = filter(res_df, !is.na(cs)),
		   aes(color = cs), shape = 21, size = 3) +
	facet_wrap(~fct_rev(stat), ncol = 1, scales = "free_y") +
	theme(axis.text.x = element_blank(),
	      axis.text.y = element_text(margin = margin(0, -.75, 0, 0, unit = "cm")),
	      axis.ticks = element_blank(),
	      panel.background = element_rect(fill = "grey96"),
	      legend.position = "top",
	      plot.margin = margin(0, 0.5, 0.5, 0.5, unit = "cm")) +
	labs(x = NULL,
	     y = NULL, 
	     color = "Credible Set:",
	     title = unique(res_df$locus))

    arrow_plot <- 
	ggplot(tracks_df |> mutate(rowid = factor(rowid))) +
	#geom_vline(data = max_pip, aes(xintercept = pos), linetype = 2, color = "grey35") +
	geom_segment(data = filter(tracks_df, feature == "intron"),
		     aes(x = start, xend = end, y = rowid, yend = rowid),
		     linewidth = .5, color = "midnightblue") +
	geom_segment(data = filter(tracks_df, feature == "exon"),
		     aes(x = start, xend = end, y = rowid, yend = rowid),
		     linewidth = 3, 
		     color = "midnightblue") +
	geom_text_repel(data = genes_df, 
		  aes(x = pos, y = rowid, label = locus),
		  fontface = "italic", 
		  size = 2.5,
		  nudge_y = .35,
		  direction = "x",
		  min.segment.length = 0,
		  segment.size = .25,
		  color = "black") +
	scale_x_continuous(limits = range(res_df$pos), labels = function(x) x/1e6L) +
	scale_y_discrete(breaks = c(0:5)) +
	theme(panel.grid = element_blank(),
	      plot.margin = margin(0, 0.5, 0.5, 0.5, unit = "cm"),
	      panel.background = element_rect(fill = "white", color = "white")) +
	labs(x = sprintf("Position in %s (Mb)", unique(res_df$chr)), 
	     y = NULL)

	ggsave("./plots/test.png", 
	       plot_grid(get_legend(pip_plot), 
			 pip_plot + theme(legend.position = "none"),
			 arrow_plot, 
			 ncol = 1, rel_heights = c(.1, 1, .35)) +
	       theme(plot.background = element_rect(fill = "white", color = "white")),
	       height = 7, width = 8)

}

## Data
# Genomics windows
regions <- 
    "./data/regions_bentham_1mb.tsv" |>
    read_tsv(col_names = c("region", "locus")) |>
    extract(region, c("chr", "left", "right"), 
	    "(chr[^:]+):(\\d+)-(\\d+)", 
	    convert = TRUE, remove = FALSE)

# GWAS data
sample_size <- 4036 + 6959
summ_stats_all <- read_tsv("./data/bentham_opengwas_1MbWindows_hg38_summstats.tsv")

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

# Analysis
#plan(multisession, workers = availableCores())
#susie_results <- future_map(regions$locus, run_susie)

susie_results <- read_rds("./susie_results.rds")

pip_df <- 
    map_dfr(setNames(susie_results, regions$locus), 
	    make_pip_df, 
	    .id = "locus")

summ_stats <- read_tsv("./data/bentham_opengwas_1MbWindows_hg38_summstats.tsv")

plot_df <- 
    left_join(pip_df, summ_stats, 
	      join_by(locus, chr, pos == bp, ref, alt)) |>
    select(locus, chr, pos, logp, pip, cs, coverage) |>
    pivot_longer(logp:pip, names_to = "stat") |>
    mutate(locus = factor(locus, levels = regions$locus),
	   stat = recode(stat, pip = "Susie PIP", logp = "GWAS p-value"),
	   stat = factor(stat, levels = c("Susie PIP", "GWAS p-value")))

plot_list <- plot_df |>
    group_split(locus) |>
    map(plot_susie)

ggsave(
   filename = "./plots/susie.pdf", 
   plot = marrangeGrob(plot_list, nrow = 1, ncol = 1), 
   width = 7, height = 5
)
    

