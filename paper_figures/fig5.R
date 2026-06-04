# ==============================================================================
# Description: Generates Figure 5 (TNFSF4 Locus and Validation).
#              Panel A: LocusZoom plot of SLE GWAS fine-mapping (SuSiE PIPs).
#              Panel B: ATAC-seq chromatin accessibility tracks at the locus.
#              Panel C: TNFSF4 low-input RNA-seq expression timecourse.
#              Panel D: CRISPRi qPCR validation of the enhancer element.
# ==============================================================================

library(tidyverse)
library(DESeq2)
library(ensembldb)
library(AnnotationHub)
library(locuszoomr)
library(rtracklayer)
library(glue)
library(ggforce)
library(ggbeeswarm)
library(cowplot)
library(ggh4x)
library(ggrepel)
library(conflicted)

ah <- AnnotationHub()
ens_data <- ah[["AH98047"]]

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::slice)

# -----------------------------------------------------------------------------
# GLOBAL AESTHETICS & SETTINGS
# -----------------------------------------------------------------------------

# Set global ggplot theme to enforce max 7pt font size
theme_set(theme_minimal(base_size = 7))

# Universal helper function to guarantee perfect title/label alignment
# Anchors text to the absolute top (y=1) and pushes it 15pt right to clear the letter label
create_title <- function(text) {
    ggdraw() + 
    draw_label(text, x = 0, y = 1, vjust = 1, hjust = 0, size = 7) + 
    theme(plot.margin = margin(t = 5, b = 7, l = 15, unit = "pt"))
}

stim_colors <- 
    read_tsv("./figure_colors_final.txt") |>
    mutate(stim = glue("{Condition} {Time}h")) |>
    select(stim, Hex) |>
    deframe()

# -----------------------------------------------------------------------------
# Fig A & B: Locus zoom + ATAC-seq
# -----------------------------------------------------------------------------
gwas_data <-
    "../gwas_finemapping/test/data/langefeld_gwas_grch38.vcf.gz" |>
    read_tsv(comment = "##") |>
    separate_rows(FORMAT:stats, sep = ":") |>
    pivot_wider(names_from = FORMAT, values_from = stats) |>
    mutate(p = 10^(-as.numeric(LP))) |>
    select(chr = 1, pos = 2, snp_id = 3, REF, ALT, p)

# In-sample LD results:
#fine_map <- 
#    "../gwas_finemapping/1-main/data/SuSiE results.SLE EA Immunochip.16OCT2024.xlsx" |>
#    readxl::read_excel() |>
#    filter(grepl("TNFSF4", region)) |>
#    select(rsid = snp, cs, pip)

rsid_df <- 
    "../gwas_finemapping/data/langefeld_munge.tsv" |> 
    read_tsv() |>
    mutate(snp_id = glue("{CHR}:{BP}:{A1}:{A2}")) |>
    select(snp_id, rsid = SNP)

fine_map <- 
    "../gwas_finemapping/results/susie_pip_TNFSF4.tsv" |>
    read_tsv() |>
    filter(cs == "L1") |>
    left_join(rsid_df)

plot_data <- 
    gwas_data |>
    left_join(fine_map, join_by(snp_id)) |>
    mutate(pip_interv = case_when(pip > 0 & pip <= 0.2 ~ "0.2",
				  pip > 0.2 & pip <= 0.4 ~ "0.4",
				  pip > 0.4 & pip <= 0.6 ~ "0.6",
				  pip > 0.6 & pip <= 0.8 ~ "0.8",
				  pip > 0.8 & pip <= 1 ~ "1"),
	   pip_interv = factor(pip_interv, levels = c("0.2", "0.4", "0.6", "0.8", "1")))

top_var <- 
    plot_data |> 
    slice_min(p, n = 1) |>
    pull(pos)

top_window <- c(top_var - 1e3, top_var + 1e3)

loc <- 
    locus(data = data.table::as.data.table(plot_data) |> mutate(chr = str_remove(chr, "chr")), 
	  gene = 'TNFSF4',
	  flank = c(0.5e5, 1e5),
	  ens_db = ens_data)

gwas_plot <- 
    ggplot(plot_data, aes(x = pos, y = -log10(p))) +
    annotate("segment",
	     x = top_var, xend = top_var, y = -Inf, yend = Inf,
	     linetype = "13", linewidth = .33) +
    geom_point(data = ~ filter(., is.na(cs)), 
	       color = "grey60") +
    geom_point(data = ~ filter(., cs == "L1"), 
	       aes(color = pip), 
	       size = 2) +
    geom_text_repel(data = ~ filter(., cs == "L1"),
		    aes(label = rsid), 
		    size = 7/.pt, box.padding = .5, 
		    nudge_x = -5000, direction = "x",
		    min.segment.length = 0) +
    scale_color_stepsn(name = "PIP:",
		       limits = c(0, 1), 
		       breaks = c(0, .2, .4, .6, .8, 1),
		       colors = c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D", "#67000D")
		       ) +
    scale_x_continuous(limits = loc$xrange,
		       labels = function(x) x/1e6L,
		       expand = c(0, 0)) +
    theme(panel.grid = element_blank(),
	  axis.title.x = element_blank(),
	  axis.text.x = element_blank(),
	  axis.title.y = element_text(size = 7, margin = margin(r = 0, unit = "lines")),
	  axis.text.y = element_text(size = 7, margin = margin(r = 0, unit = "lines")),
	  legend.text = element_text(size = 7),
	  legend.title = element_text(size = 7),
	  legend.position.inside = c(0.15, 0.7)
	  ) +
    guides(color = guide_colorbar(position = "inside", 
				  direction = "horizontal",
				  barheight = .3,
				  title.position = "top")) +
    coord_cartesian(clip = "off")

## Gene tracks
g <- 
    gg_genetracks(loc, cex.text = 7/12, 
		  filter_gene_name = loc$TX$symbol[loc$TX$symbol != ""]) +
    geom_segment(x = top_var/1e6, xend = top_var/1e6, y = -Inf, yend = Inf,
		 linetype = "13", linewidth = .33) +
    scale_x_continuous(limits = loc$xrange/1e6,
		       labels = function(x) x/1e6L,
		       expand = c(0, 0)) +
    theme(axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.line.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  panel.grid = element_blank())

## ATAC-seq
bigwigs <- 
    "../03_atacseq/1-processing/results/bwa/merged_replicate/bigwig" |>
    list.files(full.name = TRUE, pattern = "*.bigWig")

names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", basename(bigwigs))

interv <- 
    GRanges(glue("chr{loc$seqname}"), IRanges(loc$xrange[1], loc$xrange[2]))

gr <- 
    map(bigwigs, ~import(., which = interv)) |>
    map_dfr(~as.data.frame(.) |> as_tibble(), .id = "stim") |>
    mutate(stim = str_replace(stim, "_", " "),
	   stim = str_replace(stim, "unst", "Unstim"),
	   stim = paste0(stim, "h"),
	   stim = str_replace(stim, "IL4", "IL-4c"),
	   stim = str_replace(stim, "TLR7", "TLR7c"),
	   stim = str_replace(stim, "BCR", "BCRc"),
	   stim = str_replace(stim, "DN2", "DN2c"),
	   stim = factor(stim, levels = names(stim_colors)))

gr_df <- 
    gr |>
    select(stim, start, end, score) |>
    pivot_longer(c(start, end), names_to = "boundary", values_to = "bp") |>
    arrange(stim, bp)

atac_colors <- stim_colors[c("Unstim 0h", "Unstim 24h", "IL-4c 24h", "TLR7c 24h", "BCRc 24h", "DN2c 24h")]
atac_colors[c("TLR7c 24h", "BCRc 24h", "DN2c 24h")] <- stim_colors[c("TLR7c 24h", "BCRc 48h", "DN2c 48h")] 

atac_plot <- 
    ggplot(gr_df) +
    geom_ribbon(aes(x = bp, ymin = 0, ymax = score, fill = stim, color = stim), 
		linewidth = .1, outline.type = 'full', alpha = 1) +
    geom_segment(data = tibble(stim = factor("DN2c 24h", levels = levels(gr$stim))),
		 x = top_var, xend = top_var, 
		 y = max(gr$score) * 1.05, yend = -max(gr$score) * 1.9,
		 linetype = "13", linewidth = .33) +
    scale_x_continuous(limits = loc$xrange, 
		       labels = function(x) round(x/1e6L, 2),
		       expand = c(0, 0)) +
    scale_color_manual(values = c("Unstim 0h" = "black", atac_colors[-1])) +
    scale_fill_manual(values = atac_colors) +
    facet_wrap(~fct_rev(stim), ncol = 1, strip.position = "right") +
    labs(x = "Chromosome 1 (Mb)", y = "ATAC-seq peak score") +
    theme(axis.text.x = element_text(size = 7),
	  axis.text.y = element_blank(),
	  axis.title.x = element_text(size = 7),
	  axis.title.y = element_text(size = 7, hjust = .3, margin = margin(r = 0.5, unit = "lines")),
	  legend.position = "none",
	  strip.text.y.right = element_text(size = 7, angle = 0, hjust = 1, 
					    margin = margin(t = 6, l = -5, unit = "lines")),
	  strip.clip = "off",
	  panel.grid = element_blank(),
	  panel.spacing.y = unit(-3, "lines")) +
    annotate(geom = "text", x = 173.14e6, y = 0.75, size = 2.5,
	     label = glue("[0–{round(max(gr_df$score), 2)}]")) +
    coord_cartesian(clip = "off", ylim = c(0, max(gr$score)))

# -----------------------------------------------------------------------------
# Fig C: Low-input RNA-seq
# -----------------------------------------------------------------------------
gene_meta <- 
    "../01_rnaseq_lowinput/1_quantification/data/gencode.v38.primary_assembly.annotation.gtf.gz" |>
    rtracklayer::import(feature.type = "gene") |>
    as_tibble() |>
    select(gene_id, gene_name) |>
    filter(gene_name == "TNFSF4") |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

txi <- 
    read_rds("../01_rnaseq_lowinput/1_quantification/results/gene_expression_txi.rds")

tpm_df <- 
    txi$abundance |>
    as_tibble(rownames = "gene_id") |>
    inner_join(gene_meta, join_by(gene_id)) |>
    select(gene_id, gene_name, everything()) |>
    pivot_longer(-c(gene_id, gene_name), names_to = 'sample_id', values_to = 'tpm') |>
    separate(sample_id, c('donor_id', 'stim', 'timep'), sep = '_') |>
    filter(stim %in% c('Unstim', 'IL4', 'TLR7', 'BCR', 'DN2')) |>
    mutate(timep = parse_number(timep)) |>
    mutate(stim = recode(stim, 
			 'IL4' = 'IL-4c', 
			 'TLR7' = 'TLR7c',
			 'BCR' = 'BCRc', 
			 'DN2' = 'DN2c')) |>
    unite('condition', c(stim, timep), sep = ' ', remove = FALSE) |>
    mutate(stim = factor(stim, levels = c('Unstim', 'IL-4c', 'TLR7c', 'BCRc', 'DN2c')),
	   timep = factor(timep, levels = c(0, 4, 24, 48, 72)),
	   condition = paste0(condition, 'h'),
	   condition = factor(condition, levels = names(stim_colors)))

lowinp_plot <- 
    ggplot(tpm_df, aes(x = timep, y = tpm)) +
    geom_quasirandom(aes(fill = condition),
		     method = "smiley", width = .2, 
		     shape = 21, stroke = .1, size = 2) +
    scale_fill_manual(values = stim_colors) +
    facet_grid(.~stim, scale = 'free', space = 'free_x') +
    theme(axis.text = element_text(size = 7),
	  axis.title = element_text(size = 7),
	  strip.text = element_text(size = 7),
	  panel.spacing.x = unit(0.2, 'lines'),
	  panel.grid.major.y = element_blank(),
	  panel.grid.minor.y = element_blank()) +
    guides(fill = 'none') +
    labs(x = 'Time (hours)', y = 'TPM')

# -----------------------------------------------------------------------------
# Fig D: CRISPR
# -----------------------------------------------------------------------------
crispr_1 <- 
    readxl::read_excel("../supp_data/data/OX40L qPCR raw data.xlsx", 1) |>
    rownames_to_column("sample_id") |>
    pivot_longer(-sample_id, names_to = "guide")

crispr_2 <- 
    readxl::read_excel("../supp_data/data/OX40L qPCR raw data.xlsx", 2) |>
    rownames_to_column("sample_id") |>
    pivot_longer(-sample_id, names_to = "guide")
   
crispr_df <- 
    bind_rows("Cas9-KRAB" = crispr_1, "Cas9-KRAB-MeCP2" = crispr_2, .id = "version") |>
    mutate(guide = str_remove(guide, "_")) |>
    mutate_at(vars(version, guide), fct_inorder) |>
    mutate(control_test = ifelse(guide == "SCR", "control", "test"))

fig_d <- 
    ggplot(crispr_df, aes(x = guide, y = value)) +
    geom_quasirandom(aes(fill = control_test),
		     method = "smiley", width = .2, 
		     shape = 21, stroke = .5, size = 1.5) +
    geom_text(data = tribble(~version, ~x, ~y, ~lab,
			     "Cas9-KRAB", 1.5, 1.25, "****",
			     "Cas9-KRAB", 2, 1.45, "**",
			     "Cas9-KRAB-MeCP2", 1.5, 1.25, "**",
			     "Cas9-KRAB-MeCP2", 2, 1.45, "****"),
	      aes(x = x, y = y, label = lab), size = 7, size.unit = "pt") +
    geom_segment(x = 1, xend = 2, y = 1.2, yend = 1.2) +
    geom_segment(x = 1, xend = 3, y = 1.4, yend = 1.4) +
    scale_y_continuous(limits = c(0, 1.5), breaks = c(0, .5, 1)) +
    scale_fill_manual(values = c("control" = "white", "test" = "black")) +
    facet_wrap(~version, nrow = 1) +
    theme(axis.text.x = element_text(size = 7, angle = 30, hjust = 1, vjust = 1.5),
	  axis.text.y = element_text(size = 7),
	  panel.grid.minor.y = element_blank(),
	  axis.title = element_text(size = 7),
	  strip.text = element_text(size = 7),
	  strip.clip = "off",
	  plot.background = element_rect(fill = "transparent", color = NA)) +
    labs(x = NULL, 
	 y = expression(paste("  ", italic("TNFSF4"), "/", italic("GAPDH"), " expression"))) +
    guides(fill = "none")


# -----------------------------------------------------------------------------
# STANDARDIZED ASSEMBLY
# -----------------------------------------------------------------------------
fig_ab_title <- create_title("SLE risk variant at the TNFSF4/OX40L locus occurs at a DN2c-specific open chromatin region")
fig_c_title <- create_title("TNFSF4/OX40L gene expression is induced in late time points after BCR stimulation")
fig_d_title <- create_title("CRISPR inactivation near rs2205960 affects\nTNFSF4 expression")

# 1. Indent plots 
gwas_indented   <- gwas_plot + theme(plot.margin = margin(.5, .5, 0, .5, "lines"))
g_indented      <- g + theme(plot.margin = margin(0, .5, 0, .5, "lines")) 
atac_indented   <- atac_plot + theme(plot.margin = margin(0, .5, 0, .5, "lines"))
lowinp_indented <- lowinp_plot + theme(plot.margin = margin(1, 0, 0, .5, "lines"))
fig_d_indented  <- fig_d + theme(plot.margin = margin(1, .5, 0, 0, "lines"))

# 2. Assemble the A/B tracks with align="v" and axis="lr"
tracks_grid <- plot_grid(gwas_indented, g_indented, atac_indented, 
                         ncol = 1, align = "v", axis = "lr", rel_heights = c(.4, .1, 1),
                         labels = c("", "" , "b"), label_size = 10, label_y = 1, vjust = 1)

fig_ab_grid <- plot_grid(fig_ab_title, tracks_grid, ncol = 1, rel_heights = c(.05, 1),
			 labels = c("a", ""), label_size = 10, label_y = 1, vjust = 1)

# 3. Assemble C & D panels
fig_c_grid <- plot_grid(fig_c_title, lowinp_indented, ncol = 1, rel_heights = c(.1, 1),
                        labels = "c", label_size = 10, label_y = 1, vjust = 1)

fig_d_grid <- plot_grid(fig_d_title, fig_d_indented, ncol = 1, rel_heights = c(.1, 1),
                        labels = "d", label_size = 10, label_y = 1, vjust = 1)

bottom_panel <- plot_grid(fig_c_grid, NULL, fig_d_grid, nrow = 1, rel_widths = c(1, .05, .66))

# 4. Final Combination
final_figure <- 
    plot_grid(fig_ab_grid, NULL, bottom_panel, ncol = 1, rel_heights = c(1, .025, .4)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./pdf/fig5.pdf", 
       final_figure,
       width = 179,
       height = 153,
       units = "mm",
       dpi = 600,
       device = cairo_pdf)
