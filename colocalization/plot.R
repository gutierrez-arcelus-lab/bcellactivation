library(tidyverse)
library(ggsci)

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


# Colocalization
tissues <- c("macrophage", "monocyte", "neutrophil", "CD4+ T cell", "CD8+ T cell",
             "B cell", "LCL", "T cell", "blood", "Tfh cell", "Th1 cell", "Th2 cell",
             "Treg naive", "Treg memory", "CD16+ monocyte", "NK cell")

bentham_ikbke <- read_tsv("./results/bentham_region5.tsv")
bentham_ikzf1 <- read_tsv("./results/bentham_region22.tsv")
langefeld_ikbke <- read_tsv("./results/langefeld_region22.tsv")
langefeld_ikzf1 <- read_tsv("./results/langefeld_region31.tsv")

study_df <- "./data/coloc_input/eqtl_catalogue_paths.tsv" |>
    read_tsv() |>
    rename("author" = "study") |>
    mutate(study = basename(ftp_path),
	   study = sub("\\.all\\.tsv\\.gz$", "", study)) |>
    select(study, author, tissue = tissue_label, condition = condition_label, method = quant_method)

#IKZF1
ikzf1_df <- 
    bind_rows("Langefeld et al." = langefeld_ikzf1, 
	      "Bentham et al." = bentham_ikzf1, 
	      .id = "gwas") |>
    left_join(study_df, by = "study") |>
    filter(tissue %in% tissues) |>
    mutate(author = sub("_", " ", author),
	   study_label = paste0(author, " (", tissue, " ", condition, ")"),
	   method = recode(method, "tx" = "transcript", "txrev" = "splicing event"),
	   method = factor(method, levels = c("exon", "transcript", "splicing event"))) |>
    select(gwas, study = study_label, gene_id, gene_name, molecular_trait_id, method, h4)

ikzf1_df_filt <- ikzf1_df |>
    group_by(gwas, study, gene_id, method) |>
    filter(h4 == max(h4)) |>
    group_by(gwas, study, method) |>
    filter(any(h4 > .5)) |>
    group_by(gene_name) |>
    mutate(i = as.integer(any(h4 > 0.5))) |>
    ungroup() |>
    mutate(gene_name = ifelse(i == 1, gene_name, "Other"))

ikzf1_genes <- ikzf1_df_filt |> distinct(gene_name) |> pull(gene_name) |> sort()

npg_colors <- pal_npg()(10)

ikzf1_colors <- c(npg_colors[1:2], "black", "grey70")
names(ikzf1_colors) <- ikzf1_genes

ikzf1_plot <- ggplot(ikzf1_df_filt, aes(h4, study, fill = gene_name)) +
    geom_vline(xintercept = .8, linetype = 2, size = .25) +
    geom_point(size = 3, shape = 21, stroke = .2) +
    scale_fill_manual(values = ikzf1_colors) +
    scale_x_continuous(breaks = c(0, 1), limits = c(0, 1)) +
    facet_grid(method~gwas, scales = "free_y", space = "free") +
    theme(text = element_text(size = 12),
	  panel.background = element_rect(fill = "grey96")) +
    labs(x = "Probability of shared signal", 
	 y = NULL, 
	 fill = "Gene")

ggsave("./ikzf1_colocs.png", ikzf1_plot, height = 5)

#IKBKE
ikbke_df <- 
    bind_rows("Langefeld et al." = langefeld_ikbke, 
	      "Bentham et al." = bentham_ikbke, 
	      .id = "gwas") |>
    left_join(study_df, by = "study") |>
    mutate(author = sub("_", " ", author),
	   study_label = paste0(author, " (", tissue, " ", condition, ")")) |>
    select(gwas, study = study_label, gene_id, gene_name, molecular_trait_id, method, h4)

ikbke_df_filt <- ikbke_df |>
    group_by(gwas, study, gene_id, method) |>
    filter(h4 == max(h4)) |>
    group_by(gwas, study, method) |>
    filter(any(h4 > .6)) |>
    ungroup() |>
    group_by(gene_name) |>
    mutate(i = as.integer(any(h4 > 0.5))) |>
    ungroup() |>
    mutate(gene_name = ifelse(i == 1, gene_name, "Other"))

ikbke_colors <- ikbke_df_filt |> 
    distinct(gene_name) |>
    mutate(gene_name = factor(gene_name),
	   gene_name = fct_relevel(gene_name, "Other", after = Inf)) |>
    arrange(gene_name) |>
    mutate(color = c(npg_colors[c(1:2, 5)], "black", npg_colors[c(6, 9:10)], "grey80")) |>
    deframe()


ikbke_plot <- ggplot(ikbke_df_filt, aes(h4, study, fill = gene_name)) +
    geom_vline(xintercept = .8, linetype = 2, size = .25) +
    geom_point(size = 3, shape = 21, stroke = .2) +
    scale_x_continuous(breaks = c(0, 1), limits = c(0, 1)) +
    scale_fill_manual(values = ikbke_colors) +
    facet_grid(method~gwas, scales = "free_y", space = "free") +
    theme(text = element_text(size = 12),
	  panel.background = element_rect(fill = "grey96")) +
    labs(x = "Probability of shared signal", 
	 y = NULL, 
	 fill = "Gene")

ggsave("./ikbke_colocs.png", ikbke_plot)






