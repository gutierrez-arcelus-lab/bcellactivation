library(tidyverse)
library(qvalue)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

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

# fix sample mislabeling ######################################################
fix_meta <- "./mbv/fixed_ase/data/fix_metadata.tsv" |>
    read_tsv(col_names = FALSE, col_types = c(.default = "c")) |>
    extract(X4, "old_id", "\\d+_(\\d+)") |>
    select(subject_id = X1, sample_id = X2, stim = X3, old_id)

ase_fix_df <- 
    "./mbv/fixed_ase/results/ase/%s.asereadcounter.txt" |>
    sprintf(paste(fix_meta$sample_id, fix_meta$stim, sep = "_")) |>
    setNames(paste(fix_meta$sample_id, fix_meta$stim, sep = "_")) |>
    map_df(read_tsv, .id = "id") |>
    separate(id, c("sample_id", "stim"), sep = "_") |>
    extract(sample_id, c("subject_id"), "([^\\.]+)", remove = FALSE) |>
    select(subject_id, sample_id, stim, everything())

ase_df <- ase_df |>
    anti_join(fix_meta, join_by(subject_id == old_id, stim)) |>
    bind_rows(ase_fix_df)

##############################################################################

ase_res <- ase_df |>
    select(subject_id, sample_id, stim, chr = contig, position, variantID, 
	   refAllele, altAllele, refCount, total = rawDepth) |>
    mutate(chr = factor(chr, levels = paste0("chr", c(1:22, "X")))) |>
    mutate(p = pmap_dbl(list(refCount, total, 0.5), 
			function(x, y, z) binom.test(x, y, p = z, alternative = "two.sided")$p.value),
	   q = qvalue(p)$qvalues)

ase_res <- ase_res |> mutate(stim = recode(stim, "unstday0" = "unstim"))

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


annotations <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf") |>
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc")

gene_df <- annotations |>
    filter(X3 == "gene", X1 %in% paste0("chr", c(1:22, "X"))) |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(chr = X1, start = X4, end = X5, gene_id, gene_name)

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
ref_ratios <- ase_df |>
    select(sample_id, stim, var_id = variantID, refCount, totalCount, otherBases, rawDepth) |>
    mutate(total = totalCount + otherBases,
	   ref_r = refCount / total) |>
    select(sample_id, stim, var_id, ref_r) |>
    mutate(stim = factor(stim, levels = c("unstday0", "BCR", "TLR7", "DN2")))
  
ref_r_plot <- ggplot(ref_ratios, aes(ref_r)) +
    geom_histogram() +
    scale_x_continuous(breaks = c(0, .5, 1)) +
    facet_wrap(sample_id~stim, ncol = 8, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  strip.background = element_rect(fill = "grey96")) +
    labs(x = "REF allele ratio", y = NULL)

ggsave("./plots/ref_r.png", ref_r_plot, height = 12, width = 12)


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

