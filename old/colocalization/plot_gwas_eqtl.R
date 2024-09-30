library(tidyverse)
library(ggrepel)

# Harmonized data from GWAS Catalog
bentham_stats <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/bentham_GRCh38.tsv.gz" |>
    read_tsv() |>
    drop_na(hm_odds_ratio) |>
    select(chr = hm_chrom, rsid = hm_rsid, pos = hm_pos, 
	   eff_allele = hm_effect_allele, other_allele = hm_other_allele,
	   odds_ratio, beta, se = standard_error, p = p_value)

bentham_stats_hg19 <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz" |>
    read_tsv() 

bentham_stats_hg19 |> filter(rsid == "rs4917014")

ikzf1 <- read_tsv("./data/bentham_tier2_hits.tsv") |>
    filter(locus == "IKZF1") |>
    mutate(chr = sub("chr", "", chr))

ikaros_region <- bentham_stats |>
    filter(chr == 7, between(pos, ikzf1$pos - 5e5, ikzf1$pos + 5e5))

ikzf1_hg19 <- bentham_stats_hg19 |> filter(rsid == "rs4917014")

ikaros_region_hg19 <- bentham_stats_hg19 |>
    filter(chrom == 7, between(pos, ikzf1_hg19$pos - 5e5, ikzf1_hg19$pos + 5e5))

ikaros_region |> arrange(p)
ikaros_region_hg19 |> arrange(p)

ikaros_plot <- ggplot(ikaros_region, aes(pos, -log10(p))) +
    geom_point() +
    geom_label_repel(data = ikaros_region |> filter(rsid == "rs4917014"),
		     aes(label = rsid), min.segment.length = 0, direction = "y", nudge_y = .5) +
    scale_x_continuous(labels = function(x) x/1e6) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
	  panel.grid.major.x = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = "Chr7 pos (Mb)")

ggsave("./ikaros_bentham.png", ikaros_plot, width = 8, height = 4)
