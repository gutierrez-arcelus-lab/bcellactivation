library(tidyverse)
library(extrafont)

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
    read_tsv("../ase_data.tsv", col_types = "ccccccii") |>
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
    mutate(zygosity = factor(zygosity, levels = c("HET", "HOM")))

res <- 
    test_data |>
    group_by(gwas_var, stim, var_id) |>
    summarise(p = wilcox.test(imb ~ zygosity, alternative = "greater", exact = FALSE)$p.value) |>
    ungroup()

sig_imb <- filter(res, p < 0.05) |>
    distinct(var_id) |>
    pull(var_id)

# Plot
stim_colors <- 
    c("Day 0" = "#898e9f",
      "BCR" = "#003967",
      "TLR7" = "#637b31",
      "DN2" = "#a82203")

fill_colors <- stim_colors |>
    enframe("stim", "color") |>
    cross_join(tibble(zygosity = c("HOM", "HET"))) |>
    mutate(color = ifelse(zygosity == "HOM", "white", color)) |>
    unite(lab, c("stim", "zygosity"), sep = "_") |>
    deframe()

plot_data <- 
    res |>
    filter(p < 0.05) |>
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



plot_imb <- function(snp_id) {
    
    plot_data_i <- filter(plot_data, ase_var == snp_id) |>
	arrange(stim, donor_id)

    pvals_i <- res |>
	filter(var_id == unique(plot_data_i$ase_var)) |>
	mutate(p = round(p, 3),
	       p = paste("p = ", p)) |>
	mutate(stim = factor(stim, levels = names(stim_colors)))

    plot_title <- 
	sprintf("Imbalance at %s (%s)\nby genotype of %s at the %s locus.",
		unique(plot_data_i$ase_var), 
		unique(plot_data_i$ase_gene_name), 
		unique(plot_data_i$gwas_snp), 
		unique(plot_data_i$gwas_locus))

    out_plot <- 
	ggplot(plot_data_i, aes(zygosity, imb)) +
	geom_jitter(aes(color = stim, fill = lab), 
		    size = 3, shape = 21, width = .1) +
	scale_y_continuous(limits = c(NA, .5)) + 
	scale_color_manual(values = stim_colors) +
	scale_fill_manual(values = fill_colors) +
	geom_text(data = pvals_i,
		  aes(label = p),
		  size = 4, x = 0.5, y = 0.49, hjust = "inward"
		  ) +
	facet_wrap(~stim, nrow = 1) +
	theme_bw() +
	theme(
	      axis.text = element_text(size = 11),
	      axis.title = element_text(size = 11),
	      strip.text = element_text(size = 11),
	      panel.grid.minor.y = element_blank(),
	      legend.position = "none",
	      plot.title = element_text(size = 11),
	      plot.margin = margin(t = 3, b = 3., r = 2, l = 2, unit = "in"),
	      ) +
	labs(x = NULL, y = "Imbalance",
	     title = plot_title)

}


ase_plot_list <- map(sig_imb, plot_imb)
names(ase_plot_list) <- str_replace_all(sig_imb, ":", "_")


walk(names(ase_plot_list), 
     ~ggsave(sprintf("./plots/%s.pdf", .x), 
	     ase_plot_list[[.x]], 
	     width = 11,
	     height = 8.5,
	     dpi = 300))

# MIR146A region
plot_data <- 
    res |>
    distinct(gwas_var, var_id) |>
    inner_join(merged_data) |>
    extract(gwas_var, c("chr", "pos"), "([^:]+):([^:]+)", 
	    convert = TRUE) |>
    left_join(select(langefeld, -p), join_by(chr, pos)) |>
    filter(gene == "PTTG1-MIR146A") |>
    left_join(distinct(ase, var_id, annot, gene_id, gene_name)) |>
    select(gwas_locus = gene, gwas_snp = snp_id, 
	   ase_var = var_id, ase_annot = annot, 
	   ase_gene_id = gene_id, ase_gene_name = gene_name,
	   stim, donor_id, zygosity, imb) |>
    unite(lab, c("stim", "zygosity"), sep = "_", remove = FALSE) |>
    mutate(stim = factor(stim, levels = names(stim_colors)),
	   zygosity = factor(zygosity, levels = c("HOM", "HET")))

ase_plot_list_mir <- map(unique(plot_data$ase_var), plot_imb)
names(ase_plot_list_mir) <- str_replace_all(names(ase_plot_list_mir), ":", "_")

walk(names(ase_plot_list_mir), 
     ~ggsave(sprintf("./plots/mir_%s.pdf", .x), 
	     ase_plot_list_mir[[.x]], 
	     width = 11,
	     height = 8.5,
	     dpi = 300))

# P-value distribution
p_hist <- 
    ggplot(res |> mutate(stim = factor(stim, levels = names(stim_colors))), 
	   aes(x = p)) +
    geom_histogram(aes(fill = stim), bins = 20) +
    scale_fill_manual(values = stim_colors) +
    scale_x_continuous(breaks = c(0, 1)) +
    facet_wrap(~stim, ncol = 2) +
    theme_minimal() +
    theme(
	  axis.text = element_text(size = 12, family = "Arial"),
	  axis.title = element_text(size = 14, family = "Arial"),
	  legend.text = element_text(family = "Arial"),
	  legend.title = element_text(family = "Arial"),
	  strip.text = element_text(family = "Arial"),
	  panel.grid = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white")
	  ) +
    labs(x = "P-value", fill = "Stim:")


ggsave("./plots/p_hist.png", height = 4)
