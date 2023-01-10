library(tidyverse)

annotations <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v19.chr_patch_hapl_scaff.annotation.gtf" |> 
    read_tsv(comment = "#", col_types = "c-cii-c-c",
             col_names = c("chr", "feature", "start", "end", "strand", "info"))

bed <- annotations |>
    filter(chr == "chr1", feature == "gene") |>
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
           gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+")) |>
    select(chr, start, end, strand, gene_id, gene_name) |>
    filter(gene_name %in% c("IKBKE", "IL10")) |>
    mutate(x = ifelse(strand == "+", start, end),
	   xend = ifelse(strand == "+", end, start))


bentham <- read_tsv("./data/bentham_chr1.tsv")

lange <- read_tsv("./data/langefeld_chr1.tsv") |>
    mutate(logp = -log10(frequentist_add_wald_pvalue_1)) |>
    select(id = rsid, pos = position, logp)


y_lim <- c(-2.5,
	   max(c(lange$logp, bentham$logp), na.rm = TRUE))


plot_b <- ggplot(bentham, aes(pos, logp)) + 
    geom_point() +
    geom_segment(data = bed, aes(x = x, xend = xend, y = -2.5, yend = -2.5),
		 arrow = arrow(length = unit(.2, "cm"))) +
    geom_label(data = bed, aes(x = x, y = -1, label = gene_name)) +
    scale_x_continuous(labels = function(x) x/1e6L) +
    scale_y_continuous(limits = y_lim) +
    labs(x = "Position in chr1 (Mb, GRCh37)",
	 y = "-log10 (p)", title = "Bentham et al.")

plot_l <- ggplot(lange, aes(pos, logp)) + 
    geom_point() +
    geom_segment(data = bed, aes(x = x, xend = xend, y = -2.5, yend = -2.5),
		 arrow = arrow(length = unit(.2, "cm"))) +
    geom_label(data = bed, aes(x = x, y = -1, label = gene_name)) +
    scale_x_continuous(labels = function(x) x/1e6L) +
    scale_y_continuous(limits = y_lim) +
    labs(x = "Position in chr1 (Mb) (GRCh37)",
	 y = "-log10 (p)", title = "Langefeld et al.")

out <- cowplot::plot_grid(plot_b, plot_l, ncol = 1)
ggsave("./gwas_chr1.png", out)
