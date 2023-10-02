library(tidyverse)
library(qvalue)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

meta <- 
    "./array_spec.tsv" |>
    read_tsv(col_names = c("donor_id", "sample_id", "stim", "mgbid"), 
	     col_types = c(.default = "c")) |>
    unite("id", c("sample_id", "stim"), sep = "_", remove = FALSE) |>
    select(id, donor_id, sample_id, stim)
    
ase_df <- 
    sprintf("./results/%s.asereadcounter.txt", meta$id) |>
    setNames(meta$id) |>
    map_df(read_tsv, .id = "id") |>
    left_join(meta, by = "id") |>
    select(donor_id, sample_id, stim, everything()) |>
    select(-id)

ase_res <- 
    ase_df |>
    janitor::clean_names() |>
    select(donor_id, sample_id, stim, chr = contig, position, 
	   variant_id, ref_allele, alt_allele,
	   ref_count, alt_count, total_count, low_mapq_depth,
	   low_base_q_depth, raw_depth, other_bases, improper_pairs) |>
    mutate(chr = factor(chr, levels = paste0("chr", c(1:22, "X")))) |>
    mutate(p_value = pmap_dbl(list(ref_count, total_count, 0.5), 
			      function(x, y, z) binom.test(x, y, p = z, alternative = "two.sided")$p.value),
	   q_value = qvalue(p_value)$qvalues) |> 
    mutate(stim = recode(stim, "unstday0" = "Day 0"))

annotations <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v39.primary_assembly.annotation.gtf") |>
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc")

gene_df <- annotations |>
    filter(X3 == "gene", X1 %in% paste0("chr", c(1:22, "X"))) |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(chr = X1, start = X4, end = X5, gene_id, gene_name)

ase_res_annot <-
    left_join(ase_res, gene_df, join_by(chr, between(position, start, end))) |>
    select(-start, -end) |>
    group_by(donor_id, sample_id, stim, chr, position, variant_id, 
	     ref_allele, alt_allele, ref_count, alt_count, total_count,
	     low_mapq_depth, low_base_q_depth, raw_depth, other_bases, 
	     improper_pairs, p_value, q_value) |>
    summarise_at(vars(gene_id, gene_name), ~paste(., collapse = "/")) |>
    ungroup()

write_tsv(ase_res_annot, "./ase_data.tsv")



# Remove suspect samples    
ase_res_final <- ase_res_annot |>
    filter(! donor_id %in% c("10044277", "10049365", "10051708", "10068703", "10073411", "10098629"))

# Compute imbalance
imb_df <- 
    ase_res_final |>
    mutate(imb = abs(0.5 - (refCount/total))) |>
    select(gene_id, gene_name, donor_id, sample_id, stim, variantID, imb, q) |>
    mutate(stim = factor(stim, levels = c("Day 0", "BCR", "TLR7", "DN2")),
	   gene_id = str_remove(gene_id, "\\.\\d+$"))
  
# SLE genes
sle_genes <- 
    "../../../colocalization/data/coloc_input/gene_annots_bentham.tsv" |>
    read_tsv()

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


susie_res <- 
    "../../../colocalization/finemap/susie_results.tsv" |>
    read_tsv() |>
    filter(!is.na(cs))

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

selected_genes <-
    c(
      "ARID5B",
      "ATG5",
      "BANK1",
      "BLK", 
      "CD44",
      "CSK",
      "DHCR7",
      "ETS1",
      "FCGR2A",
      "LYST",
      "IFIH1",
      "IKBKE",
      "IKZF1",
      "IKZF2",
      "IKZF3",
      "IL10",
      "IL12A",
      "IRAK1",
      "IRF5",
      "IRF7",
      "IRF8",
      "ITGAM",
      "ITGAX",
      "JAZF1",
      "MECP2",
      "MIR146A",
      "MIR3142HG",
      "NCF2",
      "PTPN22",
      "PRDM1",
      "PXK",
      "RAD51B",
      "SH2B3",
      "SLC15A4",
      "SOCS1",
      "SPRED2",
      "STAT1", 
      "STAT4",
      "TASL",
      "TCF7",
      "TNFAIP3",
      "TNFSF4",
      "TNIP1",
      "TREX1",
      "TYK2",
      "UHRF1BP1",
      "UBE2L3",
      "WDFY4"
    )



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

stim_colors <- 
    c("Day 0" = "grey",
      "BCR" = "cornflowerblue",
      "TLR7" = "forestgreen",
      "DN2" = "tomato3")

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

susie_vcf <- 
    "../0-genotypes/susie_variants/allchr.mgb.vcf.gz" |>
    read_tsv(comment = "##") |>
    select(chr = `#CHROM`, pos = POS, ref = REF, alt = ALT, starts_with("2")) |>
    pivot_longer(starts_with("2"), names_to = "vcf_donor_id", values_to = "gt") |>
    separate(vcf_donor_id, c("mgb_prefix", "donor_id"), sep = "-")

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

