library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggsci)
library(ggbeeswarm)
library(patchwork)

stim_colors <- 
    c("CD40L_0hrs" = "#dbdbdb", 
     "TLR9_0hrs" = "#dbdbdb", 
     "TLR7_0hrs" = "#dbdbdb",
     "BCR_0hrs" = "#dbdbdb",
     "BCR+TLR7_0hrs" = "#dbdbdb",
     "DN2_0hrs" = "#dbdbdb",
     "CD40L_4hrs" = "#ffeba6",
     "CD40L_24hrs" = "#ffde68",
     "CD40L_48hrs" = "goldenrod3",
     "CD40L_72hrs" = "goldenrod4",
     "BCR_4hrs" = "#a1c2ed",
     "BCR_24hrs" = "#6996e3",
     "BCR_48hrs" = "#4060c8",
     "BCR_72hrs" = "#0404bf",
     "TLR7_4hrs" = "#98ab76",
     "TLR7_24hrs" = "#748f46",
     "TLR7_48hrs" = "#47632a",
     "TLR7_72hrs" = "#275024",
     "TLR9_4hrs" = "#a876d9",
     "TLR9_24hrs" = "#955bd0",
     "TLR9_48hrs" = "#803ec8",
     "TLR9_72hrs" = "#691dbf",
     "BCR+TLR7_4hrs" = "#ffa7db",
     "BCR+TLR7_24hrs" = "#ff86d0",
     "BCR+TLR7_48hrs" = "#ff61c4",
     "BCR+TLR7_72hrs" = "#ff2bb8",
     "DN2_4hrs" = "#e6907a",
     "DN2_24hrs" = "#d76b51",
     "DN2_48hrs" = "#c5432a",
     "DN2_72hrs" = "#b00000")

# SLC15A4
# Gene 
cpm_slc_df <- 
    "../bcell_lowinput/data/edger_cpm_fit.tsv" |>
    read_tsv() |>
    filter(gene_name == "SLC15A4") |>
    mutate(hours = parse_number(time),
           hours = factor(hours, levels = sort(unique(hours))),
	   stim = recode(stim, "BCR-TLR7" = "BCR+TLR7"),
           condition = paste(stim, time, sep = "_"),
	   stim = factor(stim, levels = c("CD40L", "TLR7", "TLR9", "BCR", "BCR+TLR7", "DN2"))) |>
    mutate_at(vars(cpm:fit), ~2^.)

edger_slc_results <- 
    read_tsv("../bcell_lowinput/data/edger_de_genes.tsv") |>
    mutate(stim = recode(stim, "BCR-TLR7" = "BCR+TLR7")) |>
    filter(gene_name == "SLC15A4")

p_vals_slc <- edger_slc_results |>
    inner_join(distinct(cpm_slc_df, stim, gene_id, gene_name)) |>
    select(gene_id, gene_name, stim, p = PValue, fdr = FDR) |>
    arrange(gene_name, stim) |>
    left_join(summarise(cpm_slc_df, cpm = max(cpm), .by = c(gene_id, gene_name))) |>
    mutate(p_lab = format(p, format = "e", digits = 2),
	   p_lab = paste("p =", p_lab),
	   p_lab = ifelse(fdr < 0.05, paste(p_lab, "*"), p_lab))

# Transcript
cpm_slc_tx_df <-     
    "../bcell_lowinput/data/edger_transcript_cpm_fit.tsv" |>
    read_tsv() |>
    filter(gene_name == "SLC15A4") |>
    mutate(hours = parse_number(time),
           hours = factor(hours, levels = sort(unique(hours))),
	   stim = recode(stim, "BCR-TLR7" = "BCR+TLR7"),
           condition = paste(stim, time, sep = "_"),
	   stim = factor(stim, levels = c("CD40L", "TLR7", "TLR9", "BCR", "BCR+TLR7", "DN2"))) |>
    mutate_at(vars(cpm:fit), ~2^.)

edger_slc_tx_results <- 
    read_tsv("../bcell_lowinput/data/edger_de_transcripts.tsv") |>
    mutate(stim = recode(stim, "BCR-TLR7" = "BCR+TLR7")) |>
    filter(gene_name == "SLC15A4")

p_vals_slc_tx <- edger_slc_tx_results |>
    inner_join(distinct(cpm_slc_tx_df, stim, gene_id, gene_name, tx_id)) |>
    select(gene_id, gene_name, tx_id, stim, p = PValue, fdr = FDR) |>
    arrange(gene_name, stim) |>
    left_join(summarise(cpm_slc_tx_df, cpm = max(cpm), .by = c(gene_id, gene_name, tx_id))) |>
    mutate(p_lab = format(p, format = "e", digits = 2),
	   p_lab = paste("p =", p_lab),
	   p_lab = ifelse(fdr < 0.05, paste(p_lab, "*"), p_lab))

slc_plot <-
    ggplot(data = cpm_slc_df) +
	geom_quasirandom(aes(x = hours, y = cpm, fill = condition),
			 method = "smiley", width = .2, 
			 shape = 21, stroke = .2, size = 2.5) +
	geom_text(data = p_vals_slc, 
		  aes(x = 0.5, y = cpm * 1.25, label = p_lab),
		  hjust = "inward", vjust = "inward", size = 3) +
	scale_fill_manual(values = stim_colors) + 
	facet_wrap(~fct_relevel(stim, levels(cpm_slc_df)), nrow = 2) +
	theme(panel.grid.minor = element_blank(),
	      panel.grid.major.x = element_blank(),
	      panel.background = element_rect(fill = "grey96"),
	      strip.text.x = element_text(size = 10),
	      axis.text = element_text(size = 10),
	      axis.title = element_text(size = 12),
	      legend.position = "none") +
	labs(x = NULL, y = "Normalized counts")

ggsave("./plots/slc.png", slc_plot, height = 4, width = 5)

slc_tx_plot <-
    ggplot(data = cpm_slc_tx_df) +
	geom_quasirandom(aes(x = hours, y = cpm, fill = condition),
			 method = "smiley", width = .2, 
			 shape = 21, stroke = .2, size = 2.5) +
	geom_text(data = p_vals_slc_tx, 
		  aes(x = 0.5, y = cpm * 1.25, label = p_lab),
		  hjust = "inward", vjust = "inward", size = 3) +
	scale_fill_manual(values = stim_colors) + 
	facet_grid(tx_id~fct_relevel(stim, levels(cpm_slc_tx_df)),
		   scales = "free_y") +
	theme(panel.grid.minor = element_blank(),
	      panel.grid.major.x = element_blank(),
	      panel.background = element_rect(fill = "grey96"),
	      strip.text.x = element_text(size = 10),
	      strip.text.y = element_text(size = 7),
	      axis.text = element_text(size = 10),
	      axis.title = element_text(size = 12),
	      legend.position = "none") +
	labs(x = NULL, y = "Normalized counts")

ggsave("./plots/slc_tx.png", slc_tx_plot, height = 4, width = 10)

# Gene tracks
tracks_all <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/transcript_tracks_gencodev38.tsv" |>
    read_tsv() |>
    filter(transcript_id %in% unique(cpm_slc_tx_df$tx_id))

curve_df <- tracks_all |>
    distinct(transcript_id) |>
    mutate(x0 = 128809473, x1 = 128814775)

tracks_plot <- 
    ggplot(tracks_all, aes(start, end)) +
    geom_segment(aes(x = start, xend = end, 
		     y = transcript_id, yend = transcript_id,
		     group = interaction(transcript_id, feature), 
		     linewidth = feature),
		 color = "midnightblue") +
    geom_curve(data = curve_df, 
		 aes(x = x0, xend = x1, y = transcript_id, yend = transcript_id)) +
    scale_x_continuous(labels = function(x) x/1e6L) +
    scale_linewidth_manual(values = c("intron" = 1, "exon" = 4)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  panel.border = element_blank(),
	  axis.title.y = element_blank(),
	  axis.ticks.y = element_blank(),
	  legend.position = "none") +
    labs(x = "Position on chr12 (Mb) ")

ggsave("./plots/tracks.png", tracks_plot, height = 2.5, width = 8)


# Isoform expression in Scharer et al
scharer_meta <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/Scharer/RNAseq/metadata/file_description.tsv" |>
    read_tsv() |>
    select(run, status, cell_subtype) |>
    mutate(cell_subtype = sub("^.+\\((.+)\\)$", "\\1", cell_subtype),
	   status = recode(status, "Systemic lupus erythematosus (SLE)" = "SLE"))


scharer_df <- 
    "../read_mapping/scharer/compiled_expression_humangenes.tsv" |>
    read_tsv() |>
    mutate(transcript_id = sub("\\.\\d+$", "", transcript_id)) |>
    filter(transcript_id %in% unique(cpm_slc_tx_df$tx_id)) |>
    pivot_longer(-transcript_id, names_to = "run", values_to = "tpm") |>
    left_join(scharer_meta, join_by(run)) |>
    select(run, status, cell_subtype, transcript_id, tpm) 

p <- 
    ggplot(scharer_df, aes(x = transcript_id, y = tpm, color = status)) +
    geom_boxplot() +
    geom_jitter(position = position_dodge(.7)) +
    theme_bw()

ggsave("./plots/scharer_txs.png", p)
