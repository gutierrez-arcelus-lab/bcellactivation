library(tidyverse)

langefeld <- read_tsv("./data/langefeld_hits.tsv")

genot <- 
    "./data/allchr.mgb.vcf.gz" |>
    read_tsv(comment = "##") |>
    select(-(REF:FORMAT)) |>
    pivot_longer(-(`#CHROM`:ID), names_to = "donor_id", values_to = "gt") |>
    mutate(donor_id = str_remove(donor_id, "^[^-]+-")) |>
    select(chr = `#CHROM`, pos = POS, gwas_var = ID, donor_id, gt)

windows <- 
    genot |>
    distinct(chr, pos, gwas_var) |>
    mutate(start = pos - 2e5L,
	   end = pos + 2e5L) |>
    select(chr, gwas_var, start, end)

ase <- 
    read_tsv("../ase_data_clean_annotated.tsv", col_types = "ccccccii") |>
    separate(sample_id, c("donor_id", "rep"), sep = "\\.") |>
    mutate(total = refCount + altCount,
	   imb = abs(0.5 - (refCount/total))) |>
    select(donor_id, rep, stim, var_id, annot, gene_id, gene_name, imb)

window_vars <- 
    ase |>
    distinct(var_id) |>
    separate(var_id, c("chr", "pos", "ref", "alt"), remove = FALSE, convert = TRUE) |>
    inner_join(windows, join_by(chr, between(pos, start, end))) |>
    select(gwas_var, var_id) |>
    distinct(var_id, .keep_all = TRUE)

merged_data <- 
    inner_join(ase, window_vars, join_by(var_id)) |>
    left_join(genot, join_by(donor_id, gwas_var)) |>
    mutate(zygosity = case_when(gt %in% c("0|0", "1|1") ~ "HOM",
				gt %in% c("0|1", "1|0") ~ "HET",
				TRUE ~ NA_character_)) |>
    select(donor_id, rep, stim, gwas_var, var_id, imb, zygosity) |>
    group_by(donor_id, stim, gwas_var, var_id, zygosity) |>
    summarise(imb = mean(imb)) |>
    ungroup()

test_data <- 
    merged_data |>
    group_by(gwas_var, stim, var_id, zygosity) |>
    filter(n_distinct(donor_id) >= 3) |>
    group_by(gwas_var, stim, var_id) |>
    filter(all(c("HOM", "HET") %in% zygosity)) |>
    ungroup() |>
    mutate(zygosity = factor(zygosity, levels = c("HOM", "HET")))

res <- 
    test_data |>
    group_by(gwas_var, stim, var_id) |>
    summarise(p = wilcox.test(imb ~ zygosity, alternative = "greater", exact = FALSE)$p.value) |>
    ungroup()



# test alternative greater or less
x <- 
    test_data |>
    filter(var_id == "chr10:48981904:G:T") |>
    arrange(zygosity, imb) |>
    mutate(zygosity = factor(zygosity, levels = c("HET", "HOM")))

wilcox.test(x$imb ~ x$zygosity, alternative = "two.sided", exact = FALSE)


sig_imb <- filter(res, p < 0.05)



# Plot
stim_colors <- 
    c("Day 0" = "grey",
      "BCR" = "cornflowerblue",
      "TLR7" = "#09820d",
      "DN2" = "#822808")

fill_colors <- stim_colors |>
    enframe("stim", "color") |>
    cross_join(tibble(zygosity = c("HOM", "HET"))) |>
    mutate(color = ifelse(zygosity == "HOM", "white", color)) |>
    unite(lab, c("stim", "zygosity"), sep = "_") |>
    deframe()


plot_data <- 
    sig_imb |>
    distinct(gwas_var, var_id) |>
    inner_join(merged_data) |>
    extract(gwas_var, c("chr", "pos"), "([^:]+):([^:]+)", 
	    convert = TRUE) |>
    left_join(select(langefeld, -p), join_by(chr, pos)) |>
    left_join(distinct(ase, var_id, annot, gene_id, gene_name)) |>
    select(gwas_locus = gene, gwas_snp = snp_id, 
	   ase_var = var_id, ase_annot = annot, 
	   ase_gene_id = gene_id, ase_gene_name = gene_name,
	   stim, donor_id, zygosity, imb) |>
    unite(lab, c("stim", "zygosity"), sep = "_", remove = FALSE) |>
    mutate(stim = factor(stim, levels = names(stim_colors)),
	   zygosity = factor(zygosity, levels = c("HOM", "HET")))

plot_data_i <- filter(plot_data, ase_var == unique(ase_var)[1]) |>
    arrange(stim, donor_id)

pvals_i <- filter(res, var_id == unique(plot_data_i$ase_var))

test_p <- 
    ggplot(plot_data_i, aes(zygosity, imb)) +
    geom_jitter(aes(color = stim, fill = lab), 
		size = 3, shape = 21, width = .1) +
    scale_y_continuous(limits = c(NA, .5)) + 
    scale_color_manual(values = stim_colors) +
    scale_fill_manual(values = fill_colors) +
    facet_wrap(~stim, nrow = 1) +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
	  legend.position = "none") +
    labs(x = NULL, y = "Imbalance")

ggsave("./plots/test.png", test_p, width = 6, height = 2)

