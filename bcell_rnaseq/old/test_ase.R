library(tidyverse)
library(qvalue)
library(ggrepel)
library(RColorBrewer)

meta <- "./arrayspec_ase.tsv" |>
    read_tsv(col_names = c("subject_id", "sample_id", "stim", "mgbid"), 
	     col_types = c(.default = "c")) |>
    unite("id", c("sample_id", "stim"), sep = "_", remove = FALSE) |>
    select(id, subject_id, sample_id, stim)
    
ase_df <- sprintf("./results/ase/%s.asereadcounter.txt", meta$id) |>
    setNames(meta$id) |>
    map_df(read_tsv, .id = "id") |>
    left_join(meta, by = "id") |>
    select(subject_id, sample_id, stim, everything()) |>
    select(-id)

ase_res <- ase_df |>
    select(subject_id, sample_id, stim, chr = contig, position, variantID, 
	   refAllele, altAllele, refCount, total = rawDepth) |>
    mutate(chr = factor(chr, levels = paste0("chr", c(1:22, "X")))) |>
    group_by(stim) |>
    mutate(p_null = mean(refCount/total)) |>
    ungroup() |>
    mutate(p = pmap_dbl(list(refCount, total, p_null), 
			function(x, y, z) binom.test(x, y, p = z, alternative = "two.sided")$p.value),
	   q = qvalue(p)$qvalues)

# Manhattan plot
stim_colors <- c("unstday0" = "grey",
		 "BCR" = "#55bdfb",
		 "TLR7" = "#488f31",
		 "DN2" = "#de425b")

annotations <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf") |>
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc")

gene_df <- annotations |>
    filter(X3 == "gene", X1 %in% paste0("chr", c(1:22, "X"))) |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(chr = X1, start = X4, end = X5, gene_id, gene_name)

ase_cumm <- ase_res |>
    group_by(chr) |>
    summarise(max_pos = max(position)) |>
    mutate(pos_add = lag(cumsum(max_pos), default = 0)) |>
    select(chr, pos_add)

ase_data <- ase_res |> 
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
    top_n(25, -log10(p)) |>
    ungroup() |>
    mutate(chr = factor(chr, levels = paste0("chr", c(1:22, "X")))) |>
    arrange(chr, -log10(p))

threshold_line <- ase_data |>
    filter(q < 0.01) |>
    group_by(stim) |>
    slice(which.max(p)) |>
    ungroup()
    

manh_plot <- 
    ggplot(ase_data, 
	   aes(x = pos_cumm, y = -log10(p), color = chr)) +
    geom_point(size = .5, alpha = 0.25) +
    geom_point(data = ase_top, 
	       aes(x = pos_cumm, y = -log10(p)),
	       shape = 1, size = 2, alpha = 1) +
    geom_text_repel(data = ase_top, 
		    aes(x = pos_cumm, y = -log10(p), label = gene_name),
		    inherit.aes = FALSE,
		    direction = "both", angle = 90, size = 3, nudge_y = 5, fontface = "bold", 
		    min.segment.length = 0, segment.size = .1, max.overlaps = Inf) +
    geom_hline(data = threshold_line, 
	       aes(yintercept = -log10(p)),
	       linetype = 2) +
    scale_x_continuous(label = sub("chr", "", axis_set$chr), breaks = axis_set$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ase_ylim)) +
    scale_color_manual(values = rep(c("#8aa1b4", "grey35"), unique(length(axis_set$chr)))) +
    facet_wrap(~stim, ncol = 1) +
    labs(x = NULL, y = "-log10 (p-value)") + 
    theme_minimal() +
    theme(text = element_text(size = 14),
	  legend.position = "none",
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white"))

ggsave("./plots/ase_manhattan.png", manh_plot, width = 12, height = 16)




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


#mapped_reads <- read_tsv("./star_n_uniq_mapped_reads.tsv")

# IRF8

vcf <- read_tsv("./data/allchr.mgb.vcf.gz", comment = "##")

ase_irf <- ase_res |>
    select(subject_id, sample_id, stim, chr, position, variantID, refCount, total, p, q) |>
    mutate(subject_id = as.character(subject_id),
	   eff = abs(0.5 - refCount/total),
	   ref_frac = refCount/total,
	   alt_frac = 1 - ref_frac) |>
    left_join(gene_df, join_by(chr, position >= start, position <= end)) |>
    filter(gene_name == "IRF8") |>
    select(-gene_name, -gene_id, -start, -end)

vcf_irf <- vcf |>
    filter(ID %in% ase_irf$variantID) |>
    pivot_longer(-(1:9), names_to = "subject_id", values_to = "gt") |>
    select(subject_id, position = POS, variantID = ID, gt) |>
    mutate(subject_id = sub("^\\d+_[A-Z0-9]+-(\\d+)$", "\\1", subject_id)) |>
    separate(gt, c("a1", "a2"), sep = "\\|", convert = TRUE) |>
    select(subject_id, variantID, a1, a2)
    
ase_plot_df <- left_join(ase_irf, vcf_irf, by = c("subject_id", "variantID")) |>
    select(-subject_id)

irf_df <- ase_plot_df |>
    filter(stim == first(stim)) |>
    mutate(dose = map2_int(a1, a2, `+`)) |>
    arrange(position) |>
    mutate(variantID = fct_inorder(variantID))

p <- ggplot(irf_df, aes(x = variantID, y = eff)) +
    geom_point() +
    geom_hline(yintercept = 0.5, linetype = 2) +
    facet_wrap(~sample_id, ncol = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = NULL)

ggsave("./plots/irf8.png", p, width = 8, height = 20)

p <- irf_df |>
    filter(grepl("10061340", sample_id)) |>
    mutate(coverage = case_when(total >= 10 & total < 20 ~ "10-20",
				total >= 20 & total < 40 ~ "20-40",
				total >= 40 & total < 60 ~ "40-60",
				total >= 60 & total < 80 ~ "60-80",
				total >= 80 & total < 100 ~ "80-100",
				total >= 100 ~ ">100"),
	   coverage = factor(coverage, levels = c("10-20", "20-40", "40-60", "60-80", "80-100", ">100"))) |>
    ggplot(aes(x = variantID, y = ref_frac, fill = coverage)) +
    geom_point(size = 4, shape = 21, stroke = .2) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = brewer.pal("Reds", n = 6)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = NULL, y = "REF ratio", 
	 title = "3 replicates of individual 10061340 at the IRF8 gene (coverage)")

p2 <- irf_df |>
    filter(grepl("10061340", sample_id)) |>
    mutate(coverage = case_when(total >= 10 & total < 20 ~ "10-20",
				total >= 20 & total < 40 ~ "20-40",
				total >= 40 & total < 60 ~ "40-60",
				total >= 60 & total < 80 ~ "60-80",
				total >= 80 & total < 100 ~ "80-100",
				total >= 100 ~ ">100"),
	   coverage = factor(coverage, levels = c("10-20", "20-40", "40-60", "60-80", "80-100", ">100"))) |>
    ggplot(aes(x = variantID, y = ref_frac, fill = -log10(p))) +
    geom_point(size = 4, shape = 21, stroke = .2) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_gradient(low = brewer.pal("Reds", n = 6)[1],
			high = brewer.pal("Reds", n = 6)[6]) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = NULL, y = "REF ratio", 
	 title = "3 replicates of individual 10061340 at the IRF8 gene (p-value)")

ggsave("./plots/irf8.png", cowplot::plot_grid(p, p2, ncol = 1), width = 8, height = 10)


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


gene_quants <- read_tsv("./quantifications_genes.tsv")

bach2_gene_tpm <- gene_quants |> 
    filter(sample_id %in% c("10085290.1", "10044277.1"), gene_name == "BACH2") |>
    mutate(stim = factor(stim, levels = c("unstday0", "BCR", "TLR7", "DN2"))) |>
    ggplot(aes(x = stim, y = tpm, fill = stim)) +
	geom_col() +
	scale_fill_manual(values = stim_colors) +
	facet_wrap(~sample_id, ncol = 1) +
	theme_bw() +
	theme(legend.position = "none",
	      panel.grid.minor.y = element_blank()) +
	labs(x = NULL, y = "Transcripts per Million")

ggsave("./plots/bach2_tpm.png", bach2_gene_tpm, width = 4, height = 5) 
    
tx_quants <- read_tsv("./quantifications_transcripts.tsv")

bach2_tx_quants <- tx_quants |> 
    filter(sample_id %in% c("10085290.1", "10044277.1"), gene_name == "BACH2") |>
    mutate(stim = factor(stim, levels = c("unstday0", "BCR", "TLR7", "DN2")))

bach2_tx_tpm <- bach2_tx_quants |>
    ggplot(aes(x = stim, y = tpm, fill = stim)) +
	geom_col() +
	scale_fill_manual(values = stim_colors) +
	facet_grid(tx_id~sample_id) +
	theme_bw() +
	theme(legend.position = "none",
	      panel.grid.minor.y = element_blank()) +
	labs(x = NULL, y = "Transcripts per Million")

ggsave("./plots/bach2_tx_tpm.png", bach2_tx_tpm, width = 5, height = 10) 

ase_genes |> 
    count(gene_id, gene_name, variantID, sort = TRUE) |> 
    print(n = 50)

ase_genes |> 
    distinct(stim_activ, gene_name) |>
    count(stim_activ)


# OAS3

oas3 <- ase_res |>
    filter(chr == "chr12") |>
    left_join(gene_df, join_by(chr, position >= start, position <= end)) |>
    select(-start, -end) |>
    filter(gene_name == "OAS3") |>
    mutate(eff = abs(0.5 - refCount/total)) |>
    select(sample_id, stim, position, variantID, refCount, total, eff, p) |>
    group_by(variantID) |> 
    filter(n_distinct(sample_id) > 3) |> 
    group_by(sample_id, stim) |>
    filter(n_distinct(variantID) > 3) |>
    group_by(sample_id) |>
    filter("unstday0" %in% stim) |>
    ungroup() |>
    arrange(position) |>
    mutate(variantID = fct_inorder(variantID))

oas3_filt <- oas3 |>
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


oas_p <- ggplot(oas3, aes(x = factor(i), y = refCount/total, group = variantID)) +
    geom_line(linewidth = .1) +
    geom_point(aes(color = stim, size = total)) +
    scale_color_manual(values = c("unstday0" = "grey", "BCR" = "cornflowerblue", 
				  "TLR7" = "forestgreen", "DN2" = "tomato3")) +
    scale_size(range = c(0.5, 5)) +
    facet_wrap(~as.character(sample_id), nrow = 1) +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
	  panel.grid.major.x = element_blank(),
	  legend.position = "top",
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = "variants in OAS3", y = "Allelic imbalance") +
    guides(color = guide_legend(override.aes = list(size = 6)))

tx_quant <- read_tsv("./quantifications_transcripts.tsv")

oas_quant <- filter(tx_quant, 
		    gene_name == "OAS3",
		    sample_id %in% unique(oas3$sample_id)) |>
    group_by(tx_id) |>
    filter(mean(tpm > 1) > .2) |>
    ungroup() |>
    distinct(tx_id) |>
    pull(tx_id)

tx_annot <- read_rds("../splicing/plot_data/transcripts_annot.rds")

txi_df <- tx_annot |> 
    filter(transcript_id %in% oas_quant) |>
    group_by(transcript_id) |>
    mutate(i = row_number(),
           col = case_when(feature == "intron" & i == min(i) ~ NA_character_,
                           feature == "intron" & i == max(i) ~ NA_character_,
                           feature == "exon" ~ "exon",
                           TRUE ~ "intron")) |>
    ungroup() |>
    filter(!is.na(col)) |>
    mutate(transcript_id = sub("\\.\\d+$", "", transcript_id))

oas_positions <- sort(unique(oas3$position))

plot_annot_df <- tibble(x = oas_positions, y = n_distinct(txi_df$transcript_id) + 2) |>
    rownames_to_column("i")

test <- ggplot() +
    geom_vline(xintercept = oas_positions, linewidth = .1) +
    geom_segment(data = txi_df |>
		 filter(feature == "exon"), 
		 aes(x = start, xend = end, 
		     y = transcript_id, yend = transcript_id),
		 color = "midnightblue", linewidth = 2.5) +
    geom_segment(data = txi_df |>
		 filter(feature == "intron"), 
		 aes(x = start, xend = end, 
		     y = transcript_id, yend = transcript_id),
		 color = "midnightblue", linewidth = .5) +
    geom_text_repel(data = plot_annot_df, 
		    aes(x = x, y = y, label = i), 
		    force_pull = 0,
		    nudge_y = 0.1, 
		    direction = "x", 
		    min.segment.length = 0) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = "white"),
	  axis.text.x = element_blank(),
	  axis.title = element_blank(),
	  panel.grid = element_blank(),
	  plot.margin = margin(1, 0, 0, 0, unit = "cm"))

ggsave("./plots/oas3.png", 
       cowplot::plot_grid(oas_p, test, rel_heights = c(1, .5), ncol = 1),
       width = 10, height = 5)

