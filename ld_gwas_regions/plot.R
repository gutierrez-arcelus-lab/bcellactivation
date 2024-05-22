library(tidyverse)
library(extrafont)
library(glue)
library(AnnotationHub)
library(locuszoomr)
library(patchwork)
library(rtracklayer)

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
    "./data/chr2:190970120-192970120.ld" |>
    data.table::fread() |>
    as_tibble() |>
    setNames(ld_vcf$ID) |>
    add_column(snp_id = ld_vcf$ID, .before = 1)

ld_risk_var <-
    ld_plink |>
    dplyr::filter(grepl(risk_var, snp_id)) |>
    pivot_longer(-snp_id, names_to = "var_id", values_to = "r2") |>
    dplyr::select(var_id, r2) |>
    separate(var_id, c("rsid", "ref", "alt"), sep = "-", convert = TRUE)

langefeld <- 
    "./data/langefeld_stat4.tsv" |>
    data.table::fread() |>
    dplyr::mutate(chrom = "2") |>
    dplyr::select(chrom, pos, rsid, other_allele = alleleA, effect_allele = alleleB, p)

langefeld_grch38 <- 
    "./data/langefeld_hg38.bed" |>
    read_tsv(col_names = FALSE) |>
    dplyr::select(chrom = X1, pos = X3, rsid = X4) |>
    dplyr::mutate(chrom = str_remove(chrom, "chr")) |>
    distinct()

langefeld_ld <-
    langefeld |>
    left_join(ld_risk_var, join_by(rsid, other_allele == ref, effect_allele == alt)) |>
    group_by(chrom, pos, rsid) |>
    nest() |>
    ungroup() |>
    inner_join(langefeld_grch38, join_by(chrom, rsid)) |>
    dplyr::select(chrom, pos = pos.y, rsid, data) |>
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
susie_langefeld <- 
    read_tsv("../colocalization/finemap/susie_results_langefeld.tsv")

susie_langefeld_stat4 <- 
    susie_langefeld |>
    filter(locus == "STAT4")

susie_langefeld_stat4 |>
    filter(!is.na(cs))

loc_langefeld$data |>
    as_tibble() |>
    left_join(susie_langefeld_stat4, 
	      join_by(pos))


gwas_stat4_plot <- 
    ggplot(langefeld_ld, aes(pos, -log10(p))) +
    geom_point() +
    scale_x_continuous(labels = function(x) x/1e6L) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    coord_cartesian(xlim = range(loc_langefeld$data$pos)) +
    labs(x = "Chromosome 2 (Mb)",
	 y = expression("-log"["10"]("p")))

ggsave("./plots/gwas_stat4.png", gwas_stat4_plot, width = 6, height = 3)

