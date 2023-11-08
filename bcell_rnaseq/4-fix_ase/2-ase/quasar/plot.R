library(tidyverse)
library(ggbeeswarm)


if (!file.exists("plots")) dir.create("plots")

quasar_gt <- 
    read_tsv("./quasar_genotypes.tsv", col_types = "ccddddd")

mgb_gt <- 
    read_tsv("../../0-genotypes/data/allchr.mgb.vcf.gz", comment = "##") |>
    filter(ID %in% unique(quasar_gt$snp_id)) |>
    select(snp_id = ID, matches("^\\d")) |>
    pivot_longer(-snp_id, names_to = "donor_id", values_to = "gt") |>
    mutate(donor_id = str_extract(donor_id, "(?<=-)\\d+")) |>
    separate(gt, c("a1", "a2"), sep = "\\|", convert = TRUE) |>
    mutate(gt = a1 + a2) |>
    select(donor_id, snp_id, gt)

plot_df <- 
    quasar_gt |>
    pivot_longer(g0:g2, names_to = "genot_class", values_to = "gp") |>
    left_join(mgb_gt, join_by(donor_id, snp_id)) |>
    filter(genot_class == "g0" & gt == 0 |
	   genot_class == "g1" & gt == 1 |
	   genot_class == "g2" & gt == 2)
    
genot_plot <- 
    ggplot(plot_df, aes(x = factor(gt), y = gp)) +
	geom_quasirandom(method = "smiley", width = .2, size = .5, alpha = .25) +
	scale_y_continuous(breaks = c(0, 1)) +
	theme_bw() +
	theme(panel.grid = element_blank()) +
	facet_wrap(~donor_id, nrow = 4) +
	labs(x = "Genotype", y = "QuASAR genotype probability")

ggsave("./plots/genots.png", genot_plot)

# Potential false homozygotes in MGBB
left_join(quasar_gt, mgb_gt) |>
    filter(gt != 1, g1 > .99)

left_join(quasar_gt, mgb_gt) |>
    filter(gt != 1, g1 > .99) |>
    add_count(donor_id) |>
    filter(n == max(n)) |>
    select(-n)

# Potential false heterozygotes in MGBB
# But this is confounded by monoallelic expression
left_join(quasar_gt, mgb_gt) |>
    filter(gt == 1, g1 < .1)

left_join(quasar_gt, mgb_gt) |>
    filter(gt == 1, g1 < .1) |>
    separate(snp_id, c("chr", "pos", "ref", "alt"), sep = ":", convert = TRUE) |>
    filter(chr == "chr6") |>
    mutate(tile = ntile(pos, 10)) |>
    group_by(tile) |>
    summarise(n = n()) |>
    ungroup()

