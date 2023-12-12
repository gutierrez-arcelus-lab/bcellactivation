library(tidyverse)
library(qvalue)
library(VGAM)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(ggh4x)

# plotting colors
stim_colors <- 
    c("Day 0" = "grey",
      "BCR" = "cornflowerblue",
      "TLR7" = "#09820d",
      "DN2" = "#822808")

# ASE results
ase_df <- read_tsv("./ase_data.tsv", col_types = "ffccccii") 

# Binomial test
ase_res <- 
    ase_df |>
    mutate(total = refCount + altCount) |>
    mutate(p_value = map2_dbl(refCount, total, 
			      ~binom.test(.x, .y, p = .5, alternative = "two.sided")$p.value),
	   q_value = qvalue(p_value)$qvalues) |>
    select(sample_id, stim, var_id, p_value, q_value)

## Beta-binomial
#bbinom.test <- function(x, y, p, r) {
#    
#    x0 <- pmin(x, y)
#    x1 <- pmax(x, y)
#    
#    pval1 <- pbetabinom(x0, x0 + x1, prob = p, rho = r)
#    pval2 <- (1 - pbetabinom(x1 - 1L, x0 + x1, prob = p, rho = r))
#
#    pval1 + pval2
#}
#
#bbinom_params <- 
#   test <- ase_df |>
#    filter(sample_id == first(sample_id), stim == first(stim)) |>
#    {function(x) split(x, list(x$sample_id, x$stim))}() |>
#    keep(~nrow(.) > 0) |>
#    map(~select(., refCount, altCount)) |>
#    map(as.matrix) |>
#    map(~Coef(vglm(. ~ 1, betabinomialff, trace = FALSE))) |>
#    bind_rows(.id = "id") |>
#    separate(id, c("donor_id", "rep", "stim"), sep = "\\.") |>
#    unite("sample_id", c("donor_id", "rep"), sep = ".")
#
#ase_res_bb <- ase_df |>
#    left_join(bbinom_params, join_by(sample_id, stim)) |>
#    mutate(p_value = pmap_dbl(list(refCount, altCount, mu, rho), bbinom.test))
#	   
#
#    q_value = qvalue(p_value)$qvalues) |>
#    select(sample_id, stim, var_id, p_value, q_value)
    

## Summary
#ase_res |>
#    group_by(sample_id, stim) |>
#    summarise(n_sig = sum(q_value <= 0.05)) |>
#    ungroup() |>
#    mutate(stim = fct_relevel(stim, "Day 0", after = 0)) |>
#    arrange(sample_id, stim) |>
#    write_tsv("./ase_total_fdr5.tsv")




# SLE genes
sle_genes <- 
    "../../../bcell_lowinput/data/sle_curated_genes.txt" |>
    read_lines() |>
    c("IRF1")


ase_res_annot |>
    filter(q_value < 0.05) |>
    left_join(ase_df) |>
    mutate(total = refCount + altCount,
	   imb = abs(0.5 - (refCount/total))) |>
    filter(imb > .2) |>
    group_by(var_id) |>
    filter(n_distinct(sample_id) > 5) |>
    ungroup() |>
    distinct(var_id, gene_id, gene_name) |>
    distinct(gene_name) |> pull(gene_name)


# Most extreme p-values in SLE genes
top_ase_variants <- 
    ase_res_annot |>
    filter(gene_name %in% sle_genes, q_value < 0.05) |>
    group_by(gene_id, var_id) |>
    slice_min(p_value) |>
    group_by(gene_id) |>
    top_n(10, -log10(p_value)) |>
    ungroup() |>
    distinct(gene_id, gene_name, var_id)

top_vars_df <- 
    inner_join(ase_res_annot, top_ase_variants, join_by(var_id, gene_id, gene_name)) |>
    left_join(ase_df, join_by(sample_id, stim, var_id)) |>
    mutate(stim = factor(stim, levels = c("Day 0", "BCR", "TLR7", "DN2"))) |>
    select(gene_name, sample_id, stim, var_id, refCount, altCount, p_value, q_value) |>
    arrange(gene_name, var_id, stim, p_value)

# Create data.frame for plotting
# include all conditions sequenced, 
# even though data is not observed in ASEReadCounter (no expression)

all_samples <-  
    meta |>
    select(-id) |>
    mutate(stim = recode(stim, "unstday0" = "Day 0"),
	   stim = factor(stim))

top_vars_counts <- 
    top_vars_df |>
    select(-p_value, -q_value) |>
    pivot_longer(refCount:altCount, names_to = "allele", values_to = "count") |>
    mutate(allele = str_remove(allele, "Count"),
	   allele = toupper(allele),
	   allele = factor(allele, levels = c("REF", "ALT"))) |>
    group_by(gene_name, sample_id, var_id, allele) |>
    complete(stim, fill = list(count = 0)) |>
    ungroup() |>
    arrange(gene_name, var_id, sample_id, stim, allele) |>
    inner_join(all_samples, join_by(sample_id, stim))

top_vars_pvalue <- 
    top_vars_df |>
    select(-refCount, -altCount) |>
    mutate(fdr = ifelse(q_value <= 0.05, "*", ""),
	   p_value = case_when(as.numeric(p_value) == 1L ~ "1",
			       is.na(p_value) ~ NA_character_,
			       TRUE ~ format(p_value, digits = 2, scientific = TRUE)),
	   p_value = str_replace(p_value, "e", "x10^"),
	   p_value = ifelse(is.na(p_value), p_value, paste0("p = ", p_value, fdr))) |>
    select(-q_value, -fdr)

fill_colors <- stim_colors |>
    enframe("stim", "color") |>
    cross_join(tibble(allele = c("REF", "ALT"))) |>
    mutate(color = ifelse(allele == "ALT", "white", color)) |>
    unite(lab, c("stim", "allele"), sep = "_") |>
    deframe()

# Plot data
sle_variants <- 
    top_vars_counts |>
    distinct(gene_name, var_id) |>
    arrange(gene_name) |>
    pull(var_id)

plot_ase <- function(i_variant) {

    i_df <- filter(top_vars_counts, var_id == i_variant) |>
	mutate(fill_lab = paste(stim, allele, sep = "_"))

    i_pvalues <- filter(top_vars_pvalue, var_id == i_variant)
    i_title <- sprintf("%s: %s", i_df$gene_name[1], i_variant) 

    i_margin <- 11 - (n_distinct(i_df$sample_id)/15 * 11) - .5

    ggplot(i_df) +
	geom_col(aes(y = allele, x = count, color = stim, fill = fill_lab)) +
	scale_x_continuous(expand = c(0, 0),
			   breaks = scales::pretty_breaks(2)) +
	scale_color_manual(values = stim_colors) +
	scale_fill_manual(values = fill_colors) +
	geom_text(data = filter(i_pvalues, !grepl("\\*$", p_value)),
		  aes(label = p_value), x = 0, y = 3,
		  size = 3, hjust = "inward", vjust = 0.7, color = "grey40") +
	geom_text(data = filter(i_pvalues, grepl("\\*$", p_value)),
		  aes(label = p_value), x = 0, y = 3,
		  size = 3, hjust = "inward", vjust = 0.7, color = "red", fontface = "bold") +
	#facet_grid2(sample_id~stim, scales = "free_x", drop = FALSE, independent = "x") +
	facet_grid(sample_id~stim) +
	theme_minimal() +
	theme(text = element_text(size = 12),
	      legend.position = "none",
	      strip.text.x = element_text(vjust = 3, face = "bold"),
	      strip.text.y = element_text(angle = 0, face = "bold"),
	      panel.grid.major.y = element_blank(),
	      panel.spacing = unit(1, "lines"),
	      plot.margin = margin(t = i_margin/2, b = i_margin/2, r = 1, l = 1, unit = "in"),
	      plot.background = element_rect(color = "white", fill = "white")) +
	labs(x = "Read count", y = NULL, title = i_title) +
	coord_cartesian(clip = "off")
}

ase_plot_list <- map(sle_variants, plot_ase)
names(ase_plot_list) <- str_replace_all(sle_variants, ":", "_")


#dir.create("./plots/ase_plots")
walk(names(ase_plot_list), 
     ~ggsave(sprintf("./plots/ase_plots/%s.pdf", .x), 
	     ase_plot_list[[.x]], 
	     width = 8.5,
	     height = 11,
	     dpi = 300))


















# CD44
vcf <- 
    "../0-genotypes/data/allchr.mgb.vcf.gz" |>
    read_tsv(comment = "##")

susie_res <- 
    "../../../colocalization/finemap/susie_results.tsv" |>
    read_tsv() |>
    filter(!is.na(cs))

susie_vcf <- 
    "../0-genotypes/susie_variants/allchr.mgb.vcf.gz" |>
    read_tsv(comment = "##") |>
    select(chr = `#CHROM`, pos = POS, ref = REF, alt = ALT, starts_with("2")) |>
    pivot_longer(starts_with("2"), names_to = "vcf_donor_id", values_to = "gt") |>
    separate(vcf_donor_id, c("mgb_prefix", "donor_id"), sep = "-")

cd44_ase <- vcf |>
    filter(`#CHROM` == "chr11", POS == 35208126) |>
    select(ID, matches("^\\d")) |>
    pivot_longer(-ID, names_to = "donor_id", values_to = "geno") |>
    extract(donor_id, "donor_id", ".+-(.+)") |>
    select(donor_id, geno)

susie_cd44 <- susie_res |>
    filter(locus == "CD44", pos == 35073939) |>
    select(chr, pos, ref, alt, cs) |>
    inner_join(susie_vcf)
   
imb_df <- 
    ase_df |>
    filter(var_id == "chr11:35208126:T:C") |>
    separate(sample_id, c("donor_id", "replic"), sep = "\\.") |>
    mutate(total = refCount + altCount,
	   imb = abs(0.5 - (refCount/total))) |>
    select(donor_id, replic, stim, imb) |>
    mutate(stim = factor(stim, levels = c("Day 0", "BCR", "TLR7", "DN2")))

cd44_df <- 
    left_join(cd44_ase, susie_cd44) |>
    select(donor_id, geno_reg = gt, geno_cd44 = geno) |>
    left_join(imb_df) |>
    drop_na(imb)

test_p <- 
    ggplot(cd44_df, aes(x = stim, y = imb, color = stim)) +
    geom_point(size = 4) +
    scale_color_manual(values = c("Day 0" = "grey40", "BCR" = "cornflowerblue", "TLR7" = "forestgreen", "DN2" = "tomato3")) +
    facet_wrap(~geno_reg, nrow = 1) +
    theme_bw() +
    theme(legend.position = "none")    

ggsave("./plots/cd44.png", test_p, height = 4)







  

imb_df |>
    inner_join(sle_genes, join_by(gene_id, gene_name)) |>
    filter(q < 0.05) |>
    distinct(gene_id, gene_name, variantID) |>
    add_count(gene_id) |>
    arrange(desc(n)) |>
    distinct(gene_name, n) |> print(n = Inf)

test_df <- imb_df |> 
    filter(gene_name == "IRF8") |>
    group_by(variantID) |>
    filter(any(q < 0.05)) |>
    ungroup()



susie_res |> filter(locus == "BLK")
susie_res |> filter(grepl("IRF5", locus), grepl("L2", cs))

susie_res |> filter(locus == "IRF8")

test_reg <- 
    filter(susie_vcf, chr == "chr16", pos == 85985027) |>
    select(donor_id, gt) |>
    mutate(gt = recode(gt, "0|0" = "0", "1|1" = "0", "1|0" = "1", "0|1" = "1"))

test_p2 <- 
    inner_join(test_df, test_reg, join_by(donor_id)) |>
    ggplot(aes(x = variantID, y = imb, color = gt)) +
	geom_boxplot() +
	geom_point(position = position_dodge(width = .9)) +
	facet_wrap(~stim, ncol = 1) +
	theme_bw() +
	theme(panel.grid = element_blank(),
	      axis.text.x = element_text(angle = 90),
	      legend.position = "top")

ggsave("./plots/test2.png", test_p2)

null_vars <- imb_df |>
    filter(stim == "Day 0") |>
    group_by(variantID) |>
    filter(all(q > 0.1)) |>
    ungroup() |>
    distinct(variantID) |>
    pull(variantID)
   
imb_df |>
    filter(variantID %in% null_vars) |>
    group_by(variantID) |>
    filter(any(imb > .2 & q < 0.05)) |>
    ungroup() |>
    filter(gene_name %in% selected_genes) |>
    group_by(gene_name, variantID) |>
    summarise(n = sum(q < 0.05)) |>
    ungroup() |> arrange(desc(n))

imb_df |>
    filter(gene_name == "ITGAX", variantID == "chr16:31363214:C:G") |>
    print(n = Inf)
    
ase_df |>
    filter(variantID == "chr16:31363214:C:G",
	   donor_id %in% c("10048130", "10050385")) |>
    print(width = Inf)

test_p <- 
    ase_res |>
    filter(variantID == "chr16:31363214:C:G",
	   donor_id %in% c("10048130", "10050385")) |>
    mutate(alt_ref_ratio = (total - refCount)/refCount,
	   stim = factor(stim, levels = c("Day 0", "BCR", "TLR7"))) |>
    select(donor_id, stim, alt_ref_ratio, p) |>
    mutate(p = format(p, scientific = TRUE, digits = 2),
	   p = sub("\\.0e\\+00$", "", p),
	   p = paste("P =", p)) |>
    ggplot(aes(x = stim, y = alt_ref_ratio, fill = stim)) +
	geom_col() +
	geom_text(aes(label = p), size = 2, vjust = 0) +
	scale_fill_manual(values = c("Day 0" = "grey40", "BCR" = "cornflowerblue", "TLR7" = "forestgreen")) +
	facet_wrap(~donor_id, nrow = 1) +
	theme_bw() +
	theme(text = element_text(size = 10),
	      plot.title = element_text(size = 10),
	      legend.position = "none",
	      panel.grid.major.x = element_blank(),
	      panel.grid.minor.y = element_blank()) +
	labs(x = "Stimulation", y = "ALT / REF count ratio", 
	     title = "chr16:31363214:C:G at ITGAX")
   
ggsave("./plots/test.png", test_p, height = 2.5, width = 4)

















ase_res |> 
    filter(q < 0.05) |>
    distinct(stim, variantID) |>
    count(stim)

ase_res_annot |> 
    filter(gene_name == "BLK") |>
    group_by(sample_id, stim) |>
    filter(any(q < 0.05)) |>
    ungroup()

test_df <- 
    ase_res_annot |> 
    filter(gene_name == "BLK") |>
    group_by(sample_id) |> 
    filter(any(q < 0.05)) |>
    ungroup() |>
    group_by(position) |>
    nest() |>
    ungroup() |>
    rowid_to_column("i") |>
    unnest(cols = data)

test_plot <- 
    ggplot(test_df, aes(x = factor(i), y = refCount/total)) +
    geom_line(aes(group = variantID), linewidth = .1) +
    geom_point(aes(color = stim, size = total)) +
    scale_color_manual("Stim:", values = stim_colors) +
    scale_size("N reads:", range = c(0.5, 5)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, .5, 1), limits = c(0, 1)) +
    facet_wrap(~as.character(sample_id), nrow = 5, scales = "free_x") +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
	  panel.grid.major.x = element_blank(),
	  legend.position = "top",
	  plot.background = element_rect(fill = "white", color = "white")) +
    guides(color = guide_legend(override.aes = list(size = 4)))

ggsave("./plots/test2.png", test_plot, height = 10, width = 12)




# Susie fine-mapped variants
susie_res <- 
    "../../../colocalization/finemap/susie_results.tsv" |>
    read_tsv() |>
    filter(!is.na(cs), pip >= 0.1)


# Salmon gene expression
meta_data <- "../1-mapping/metadata.tsv" |>
    read_tsv(col_names = c("donor_id", "sample_id", "stim", "fq1", "fq2"),
	     col_types = c(.default = "c")) |>
    select(new_donor_id = donor_id, new_sample_id = sample_id, stim, fq1) |>
    extract(fq1, c("donor_id"), "^\\d+_([^_]+)_.+") |>
    mutate(rep_id = sub("^\\d+\\.(\\d)$", "\\1", new_sample_id),
	   sample_id = paste(donor_id, rep_id, sep = ".")) |>
    select(donor_id, sample_id, stim, new_donor_id, new_sample_id)

gene_expression <- 
    read_tsv("../../results/salmon_genes.tsv", col_types = "ccccdd") |>
    left_join(meta_data, join_by(sample_id, stim)) |>
    select(donor_id = new_donor_id, sample_id = new_sample_id,
	   stim, gene_id, gene_name, counts, tpm)


# TYK2
tyk2_expression_plot <- 
    gene_expression |>
    filter(gene_name == "TYK2") |>
    mutate(stim = recode(stim, "unstday0" = "Day 0"),
	   stim = factor(stim, levels = c("Day 0", "BCR", "TLR7", "DN2"))) |>
    ggplot(aes(x = stim, y = tpm)) +
	geom_col(aes(fill = stim)) +
	facet_wrap(~sample_id, nrow = 4) +
	scale_fill_manual(values = stim_colors) +
	theme_bw() +
	theme(axis.text.x = element_blank(),
	      axis.ticks.x = element_blank(),
	      panel.grid.minor.y = element_blank(),
	      panel.grid.major.x = element_blank(),
	      legend.position = "top",
	      plot.background = element_rect(fill = "white", color = "white")) +
    labs(y = "Transcripts per Million", x = NULL)

ggsave("./plots/tyk2_expression.png", tyk2_expression_plot, height = 5)

# BLK
blk_expression_plot <- 
    gene_expression |>
    filter(gene_name == "BLK") |>
    mutate(stim = recode(stim, "unstday0" = "Day 0"),
	   stim = factor(stim, levels = c("Day 0", "BCR", "TLR7", "DN2"))) |>
    ggplot(aes(x = stim, y = tpm)) +
	geom_col(aes(fill = stim)) +
	facet_wrap(~sample_id, nrow = 4) +
	scale_fill_manual(values = stim_colors) +
	theme_bw() +
	theme(axis.text.x = element_blank(),
	      axis.ticks.x = element_blank(),
	      panel.grid.minor.y = element_blank(),
	      panel.grid.major.x = element_blank(),
	      legend.position = "top",
	      plot.background = element_rect(fill = "white", color = "white")) +
    labs(y = "Transcripts per Million", x = NULL)

ggsave("./plots/blk_expression.png", blk_expression_plot, height = 5)

susie_blk <- filter(susie_res, locus == "BLK") |>
    select(locus, chr, pos, ref, alt, pip)

blk_vcf <- 
    susie_vcf |>
    inner_join(slice(susie_blk, 1), join_by(chr, pos, ref, alt)) |>
    select(donor_id, gt) |>
    mutate(is_het = case_when(gt %in% c("0|0", "1|1") ~ 0,
			      gt %in% c("0|1", "1|0") ~ 1,
			      TRUE ~ NA)) |>
    select(donor_id, is_het) |>
    mutate(is_het = recode(is_het, "0" = "No", "1" = "Yes"))
    

test_3 <- imb_df |>
    filter(variantID == "chr8:11494547:A:G") |>
    inner_join(blk_vcf) |>
    select(donor_id, sample_id, stim, variantID, imb, q, is_het) |>
    ggplot(aes(x = stim, y = imb)) +
	ggbeeswarm::geom_quasirandom(aes(color = stim, shape = factor(is_het)),
				     method = "smiley", width = .2, size = 3) +
	scale_color_manual(values = c("Day 0" = "grey40", "BCR" = "cornflowerblue", "TLR7" = "forestgreen", "DN2" = "tomato3")) +
	scale_shape_manual("Heterozygote\nat rs2736332:", values = c("No" = 1, "Yes" = 19)) +
	theme_bw() +
	theme(text = element_text(size = 10),
	      plot.title = element_text(size = 10)) +
	labs(x = "Stimulation", y = "Allelic imbalance", 
	     title = "Allelic imbalance at chr8:11494547:A:G (BLK)\nas a function of the heterozygosity at rs2736332") +
	guides(color = "none")
   
ggsave("./plots/test3.png", test_3, height = 3, width = 6)



#IRF 8
susie_res |>
    filter(locus == "IRF8")

susie_irf8 <- filter(susie_res, locus == "IRF8") |>
    select(locus, chr, pos, ref, alt, pip)

irf8_genos <- susie_vcf |>
    #inner_join(susie_irf8, join_by(chr, pos, ref, alt)) |>
    filter(chr == "chr16", pos == 85985027) |>
    select(donor_id, pos, gt) |>
    mutate(is_het = case_when(gt %in% c("0|0", "1|1") ~ 0L,
			      gt %in% c("1|0", "0|1") ~ 1L,
			      TRUE ~ NA_integer_))

irf8_df <- 
    ase_res_annot |> 
    filter(gene_name == "IRF8") |>
    group_by(sample_id) |> 
    filter(any(q < 0.05)) |>
    ungroup() |>
    group_by(position) |>
    nest() |>
    ungroup() |>
    rowid_to_column("i") |>
    unnest(cols = data) |>
    group_by(variantID) |>
    filter(any(q < 0.05)) |>
    ungroup()

irf8_plot <- 
    ggplot(irf8_df, aes(x = factor(i), y = refCount/total)) +
    geom_line(aes(group = variantID), linewidth = .1) +
    geom_point(aes(color = stim, size = total)) +
    scale_color_manual("Stim:", values = stim_colors) +
    scale_size("N reads:", range = c(0.5, 4)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, .5, 1), limits = c(0, 1)) +
    facet_wrap(~as.character(sample_id), ncol = 2) +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
	  panel.grid.major.x = element_blank(),
	  legend.position = "top",
	  plot.background = element_rect(fill = "white", color = "white")) +
    guides(color = guide_legend(override.aes = list(size = 4)))

ggsave("./plots/irf8.png", irf8_plot, height = 6, width = 6)


irf8_df |> filter(stim == "BCR") |>
    mutate(imb = abs(0.5 - (refCount/total))) |>
    select(donor_id, sample_id, stim, position, variantID, imb) |>
    left_join(irf8_genos, join_by(donor_id)) |>
    group_by(variantID, stim) |>
    filter(n_distinct(gt_state) > 1) |>
    group_by(variantID, stim, gt_state) |>
    summarise(mean_imb = mean(imb)) |>
    ungroup() |>
    pivot_wider(names_from = gt_state, values_from = mean_imb) |>
    mutate(d = het - hom)


ase_res_annot |>
    filter(chr == "chr16", between(position, 85985027 - 250000, 85985027 + 250000)) |>
    mutate(imb = abs(0.5 - (refCount/total))) |>
    select(donor_id, sample_id, stim, variantID, imb) |>
    left_join(irf8_genos, join_by(donor_id)) |>
    select(sample_id, stim, variantID, imb, is_het) |>
    group_by(stim, variantID) |>
    filter(n_distinct(sample_id) >= 6, 
	   min(table(is_het)) >= 3,
	   sd(is_het) > 0) |>
    summarise(r = cor(imb, is_het)) |>
    ungroup() |>
    arrange(desc(abs(r)))

ase_res_annot |>
    filter(variantID == "chr16:85779625:A:C", stim == "BCR") |>
    mutate(imb = abs(0.5 - (refCount/total))) |>
    select(donor_id, sample_id, stim, variantID, gene_name, imb) |>
    left_join(irf8_genos, join_by(donor_id)) |>
    select(sample_id, stim, variantID, gene_name, imb, is_het)














    
ase_filt <- ase_res |>
    group_by(variantID) |>
    filter(!any(stim == "unstim" & q < 0.1)) |>
    ungroup() |>
    filter(stim != "unstim")

# Manhattan plot
stim_colors <- c("unstim" = "#dbdbdb",
		 "BCR" = "#6996e3",
		 "TLR7" = "#748f46",
		 "DN2" = "#d76b51")

sle_genes <- read_tsv("../bcell_scrna/reported_genes.tsv")




ase_cumm <- ase_filt |>
    group_by(chr) |>
    summarise(max_pos = max(position)) |>
    mutate(pos_add = lag(cumsum(max_pos), default = 0)) |>
    select(chr, pos_add)

ase_data <- ase_filt |> 
    inner_join(ase_cumm) |> 
    mutate(pos_cumm = position + pos_add)

axis_set <- ase_data |> 
  group_by(chr) |>
  summarize(center = mean(pos_cumm))

ase_ylim <- ase_data |>
  filter(p == min(p)) |>
  mutate(ylim = abs(floor(log10(p))) + 5) |>
  pull(ylim)
   
ase_top <- ase_data |>
    group_by(variantID) |>
    filter(n_distinct(subject_id) > 1) |>
    group_by(stim) |>
    top_n(1000, -log10(p)) |>
    group_by(stim, variantID) |>
    slice_max(-log10(p)) |>
    ungroup() |>
    select(stim, chr, position, pos_cumm, p) |>
    mutate(chr = as.character(chr)) |>
    left_join(gene_df, join_by(chr, position >= start, position <= end)) |>
    group_by(stim, gene_id) |>
    slice_max(-log10(p)) |>
    group_by(stim, position) |>
    slice_max(-log10(p)) |>
    group_by(stim) |>
    top_n(50, -log10(p)) |>
    ungroup() |>
    mutate(chr = factor(chr, levels = paste0("chr", c(1:22, "X")))) |>
    arrange(chr, -log10(p))

ase_sle <- ase_data |>
    select(subject_id, sample_id, stim, chr, position, pos_cumm, variantID, p, q) |>
    left_join(gene_df, join_by(chr, position >= start, position <= end)) |>
    filter(gene_name %in% sle_genes$gene, q < 0.05) |>
    group_by(variantID) |>
    filter(n_distinct(subject_id) > 1) |>
    group_by(stim, variantID) |>
    slice_max(-log10(p)) |>
    group_by(stim, gene_id) |>
    slice_max(-log10(p)) |>
    ungroup()

threshold_line <- ase_data |>
    filter(q < 0.05) |>
    group_by(stim) |>
    slice(which.max(p)) |>
    ungroup()

x_labels <- sub("chr", "", axis_set$chr)
    
manh_plot <- 
    ggplot(ase_data, 
	   aes(x = pos_cumm, y = -log10(p), color = stim)) +
    geom_point(size = 2, alpha = 1) +
    geom_hline(data = threshold_line, 
	       aes(yintercept = -log10(p)),
	       linetype = 2, linewidth = 1.5) +
    scale_x_continuous(labels = x_labels,
		       breaks = axis_set$center,
		       expand = c(0.01, 0.01),
		       guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ase_ylim)) +
    scale_color_manual(values = stim_colors) +
    facet_wrap(~fct_relevel(stim, c("BCR", "TLR7", "DN2")), 
	       ncol = 1) +
    labs(x = NULL, y = "-log10 (p-value)") + 
    theme_bw() +
    theme(axis.text = element_text(size = 22),
	  axis.text.y = element_text(size = 22),
	  axis.title = element_text(size = 22),
	  strip.text = element_text(size = 22),
	  legend.position = "none",
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank())

ggsave("./plots/ase_manhattan.png", manh_plot, width = 12, height = 8, dpi = 600)




###############################################################################
# Some genes have a lot of SNPs with monoallelic expression of the REF allele
# mapping bias?
ase_genes <- ase_data |>
    filter(q < 0.05) |>
    select(sample_id, stim, chr, position, refCount, total, p) |>
    mutate(chr = as.character(chr)) |>
    left_join(gene_df, join_by(chr, position >= start, position <= end)) |>
    filter(!is.na(gene_name)) |>
    group_by(sample_id, stim, chr, position, p, refCount, total) |>
    summarise(gene = paste(gene_name, collapse = ",")) |>
    ungroup()

ase_genes |>
    count(stim, gene, sample_id, sort = TRUE) |>
    filter(n > 10) |> print(n = 30)


ase_genes |>
    mutate(imb = refCount/total)

ase_genes |>
    filter(chr == "chr1", position == 2508718)



ase_x <- ase_genes |>
    filter(stim == "TLR7", gene == "HLA-DPB1", sample_id == "10044277.1") |>
    arrange(position)

ase_genes |>
    filter(stim == "unstday0", gene == "RIPOR2", sample_id == "10085290.1") |>
    arrange(position) |>
    print(n = Inf)

vcf <- read_tsv("./data/allchr.r2filtered.mgb.vcf.gz", comment = "##")

    


###############################################################################


# SLE genes
sle_genes <- read_tsv("../bcell_scrna/reported_genes.tsv")

ase_select_df <- ase_res |>
    filter(total >= 20) |>
    group_by(sample_id, variantID) |>
    filter("unstday0" %in% stim) |>
    filter(n_distinct(stim) > 2) |>
    ungroup() |>
    mutate(eff = 0.5 - refCount/total) |>
    select(sample_id, stim, chr, position, variantID, eff)

ase_genes <- 
    left_join(filter(ase_select_df, stim != "unstday0"),
	      filter(ase_select_df, stim == "unstday0"), 
	      by = c("sample_id", "variantID", "chr", "position"), 
	      suffix = c("_activ", "_rest")) |>
    filter(between(eff_rest, -.2, .2)) |>
    mutate(delta = map2_dbl(eff_rest, eff_activ, function(x, y) diff(c(x, y)))) |>
    filter(abs(delta) > .2) |>
    group_by(variantID) |>
    filter(n_distinct(sample_id) > 1) |>
    ungroup() |>
    left_join(gene_df, join_by(chr, position >= start, position <= end)) |>
    select(-start, -end) |>
    filter(!is.na(gene_name)) |> 
    group_by(gene_id) |>
    filter(n_distinct(variantID) > 1) |>
    ungroup()

ase_genes |> filter(gene_name %in% sle_genes$gene)


bach2 <- ase_res |> 
    filter(variantID == "chr6:90206626:T:C",
	   sample_id %in% c("10085290.1", "10044277.1")) |>
    mutate_at(vars(1:2), as.character) |>
    select(sample_id, stim, variantID, refCount, total) |>
    mutate(imb = refCount/total) |>
    mutate(stim = fct_inorder(stim)) |>
    mutate(lab = paste(refCount, total, sep = "/"))

bach2_plot <- ggplot(bach2, aes(x = stim, y = imb, fill = stim)) +
    geom_col() +
    geom_text(aes(label = lab), vjust = 1) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~sample_id, ncol = 1) +
    theme_bw() +
    theme(legend.position = "none",
	  panel.grid.minor.y = element_blank()) +
    labs(x = NULL, y = "Ref allele / total")

ggsave("./plots/bach2_ase.png", bach2_plot, width = 4, height = 5) 



bach2_vcf <- vcf |> 
    filter(ID == "chr6:90206626:T:C") |>
    pivot_longer(-(1:9), names_to = "donor_id") |>
    select(donor_id, REF, ALT, value) |>
    extract(donor_id, "donor_id", "\\d+_[^-]+-(\\d+)")

bach2 |>
    mutate(donor_id = sub("\\.1$", "", sample_id)) |>
    select(donor_id, stim, refCount, total, lab) |>
    left_join(bach2_vcf, by = "donor_id")




# OAS3
oas3 <- ase_res |>
    filter(chr == "chr12") |>
    left_join(gene_df, join_by(chr, position >= start, position <= end)) |>
    select(-start, -end) |>
    filter(gene_name == "OAS3") |>
    mutate(eff = abs(0.5 - refCount/total)) |>
    select(subject_id, sample_id, stim, position, variantID, refCount, total, eff, p, q) |>
    group_by(variantID) |> 
    filter(n_distinct(subject_id) > 1) |> 
    group_by(sample_id, stim) |>
    filter(n_distinct(variantID) > 3) |>
    group_by(sample_id) |>
    filter("unstim" %in% stim) |>
    ungroup() |>
    arrange(position) |>
    mutate(variantID = fct_inorder(variantID))

oas3_filt <- oas3 |>
    group_by(sample_id) |>
    filter(n_distinct(stim) > 3 & sum(q < 0.05) > 5) |>
    group_by(sample_id, variantID) |>
    summarise(delta = diff(range(eff))) |>
    group_by(sample_id) |>
    filter(sum(delta > .2) > 3) |>
    ungroup() |>
    distinct(sample_id)

oas3 <- oas3 |>
    filter(sample_id %in% oas3_filt$sample_id) |>
    left_join(tibble(variantID = levels(oas3$variantID)) |>
	      rowid_to_column("i"))


oas_p <- 
    ggplot(oas3, aes(x = factor(i), y = refCount/total, group = variantID)) +
    geom_line(linewidth = .1) +
    geom_point(aes(color = stim, size = total)) +
    scale_color_manual(values = stim_colors) +
    scale_size(range = c(0.5, 5)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, .5, 1), limits = c(0, 1)) +
    facet_wrap(~as.character(sample_id), nrow = 1, scales = "free_x") +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
	  panel.grid.major.x = element_blank(),
	  legend.position = "top",
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = "variants in OAS3", y = "Reference allele ratio", 
	 color = "Stim:", size = "# of reads:") +
    guides(color = guide_legend(override.aes = list(size = 6)))

gene_tracks <- 
    read_tsv("../../genecode_V38_tracks.tsv") |>
    filter(gene_name == "OAS3") |>
    mutate(bp = map2(start, end, ~.x:.y)) |>
    select(feature, gene_name, i, bp) |>
    unnest(bp)

intron_tracks <- gene_tracks |>
    filter(feature == "intron") |>
    group_by(i) |>
    slice(1:100) |>
    ungroup()

oas_tmp_df <- filter(gene_tracks, feature == "exon") |>
    bind_rows(intron_tracks) |>
    arrange(bp) |>
    mutate(pos = seq_len(n()))

oas_vars <- distinct(oas3, var_idx = i, position) |>
    left_join(oas_tmp_df, join_by(position == bp))

oas_track_df <- oas_tmp_df |>
    group_by(feature, i) |>
    summarise(start = min(pos), end = max(pos)) |>
    ungroup() |>
    arrange(i)

oas_track <- 
    ggplot(oas_track_df) +
    geom_segment(aes(x = start, xend = end, y = 1, yend = 1, 
		     linewidth = feature)) +
    geom_segment(data = oas_vars |> filter(as.logical(row_number() %% 2)),
		 aes(x = pos, xend = pos, y = 1, yend = 1.25),
		 linewidth = .5) +
    geom_text_repel(data = oas_vars |> 
		    filter(as.logical(row_number() %% 2)),
		    aes(x = pos, y = 1.25, label = var_idx),
		    nudge_y = 0.1, 
		    direction = "x", 
		    min.segment.length = 0,
		    size = 3,
		    segment.size = .25) +
    geom_segment(data = oas_vars |> filter(!as.logical(row_number() %% 2)),
		 aes(x = pos, xend = pos, y = 1, yend = 1.5),
		 linewidth = .5) +
    geom_text_repel(data = oas_vars |> 
		    filter(!as.logical(row_number() %% 2)),
		    aes(x = pos, y = 1.5, label = var_idx),
		    nudge_y = 0.1, 
		    direction = "x", 
		    min.segment.length = 0,
		    size = 3,
		    segment.size = .25) +
    scale_linewidth_manual(values = c("exon" = 4, "intron" = 1)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(plot.background = element_rect(color = "white", fill = "white"),
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  panel.grid = element_blank(),
	  legend.position = "none") +
    coord_cartesian(clip = "off")

ggsave("./plots/oas3.png", 
       oas_p + oas_track + plot_layout(ncol = 1, heights = c(1, .3)),
       width = 8, height = 3.5, dpi = 600)



# check OAS3 genotypes
vcf_oas <- vcf |>
    filter(ID %in% unique(oas3$variantID)) |>
    pivot_longer(-(1:9), names_to = "donor_id") |>
    select(donor_id, ID, REF, ALT, value) |>
    extract(donor_id, "donor_id", "\\d+_[^-]+-(\\d+)")
    
oas3 |>
    mutate(donor_id = sub("\\.\\d$", "", sample_id)) |>
    select(donor_id, sample_id, stim, ID = variantID, refCount, total, i) |>
    left_join(vcf_oas, join_by(donor_id, ID)) |>
    arrange(donor_id, sample_id, i, stim) |>
    filter(i == 7) |>
    mutate(imb = refCount/total) |>
    print(n = Inf)


vcf_indel <- read_tsv("./indels/data/allchr.r2filtered.mgb.vcf.gz", comment = "##")

vcf_indel |> 
    filter(`#CHROM` == "chr12", between(POS, 112962000, 112965000)) |>
    pivot_longer(-(1:9), names_to = "sample_id") |>
    extract(sample_id, "donor_id", "\\d+_[^-]+-(\\d+)")

ase_df |>
    filter(variantID == "chr12:112963394:G:A") |>
    select(sample_id, stim, variantID, refCount, altCount, totalCount, otherBases) |>
    left_join(ase_res |>
	      filter(variantID == "chr12:112963394:G:A") |>
	      select(sample_id, stim, variantID, p)) |>
    mutate(p = round(p, 2)) |>
    print(n = Inf)

# QC plots

# Test ASE in GWAS regions
regions_df <- "../colocalization/finemap/data/regions_bentham_1mb.tsv" |>
    read_tsv(col_names = c("region", "locus")) |>
    extract(region, c("chr", "start", "end"), "(chr[^:]+):(\\d+)-(\\d+)", convert = TRUE)

ase_gwas <- ase_res |>
    mutate(chr = sub("^(chr[^:]+).+$", "\\1", variantID)) |>
    inner_join(regions_df, join_by(chr, between(position, start, end))) |>
    select(sample_id, stim, chr, position, variantID, refCount, total, p, q, locus) |>
    inner_join(gene_df, join_by(chr, between(position, start, end))) |>
    select(sample_id, stim, chr, position, variantID, refCount, total, p, q, locus, gene_id, gene_name)

write_tsv(ase_gwas, "./ase_gwasregions.tsv")

