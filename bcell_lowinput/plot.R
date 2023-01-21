library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggsci)
library(ggbeeswarm)

stims_order <- c("IL4_4hrs", "IL4_24hrs",
      paste("CD40L", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("BCR", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("TLR7", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("BCR-TLR7", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("TLR9", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("DN2", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      "NA_NA")
      
stim_colors <- c(brewer.pal(n = 3, "Greys")[-1],
		 brewer.pal(n = 9, "YlOrRd")[c(1, 3, 4, 5)],
		 brewer.pal(n = 9, "Blues")[c(2, 4, 6, 8)],
		 brewer.pal(n = 9, "Greens")[c(2, 4, 6, 8)],
		 grep("cyan", colors(), value = TRUE)[c(8, 2, 4, 6)],
		 grep("pink", colors(), value = TRUE)[c(16, 7, 1, 4)],
		 paste0("tomato", c("", 2:4)),
		 "white")

names(stim_colors) <- stims_order

cpm_df <- read_tsv("./data/edger_transcript_cpm_fit.tsv")
edger_results <- read_tsv("./data/edger_de_transcripts.tsv")

edger_results %>%
    group_by(stim, gene_id) %>%
    filter(any(FDR < 0.05)) %>%
    ungroup() %>%
    distinct(stim, gene_id) %>%
    count(stim)


select_genes <- edger_results %>%
    filter(stim == "DN2", FDR < 0.05) %>%
    filter(gene_name %in% c("NSF", "RPS24", "EIF4A2", "IRF5")) %>%
    select(stim, gene_id, gene_name, tx_id)

cpm_plot_df <- inner_join(select_genes, cpm_df) %>%
    mutate(hours = parse_number(time),
           hours = factor(hours, levels = sort(unique(hours))),
           condition = paste(stim, time, sep = "_")) %>%
    mutate_at(vars(cpm:fit), ~2^.) %>%
    unite("lab", c(gene_name, tx_id), sep = "\n")


tcourse <- cpm_plot_df %>%
    ggplot(aes(x = hours, y = cpm, fill = condition)) +
    geom_quasirandom(method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 3) +
    scale_fill_manual(values = stim_colors) + 
    facet_wrap(~lab, scales = "free_y") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
	  panel.grid.major.x = element_blank(),
	  legend.position = "none") +
    labs(x = NULL)

ggsave("./plots/timecourse_tx.png", tcourse) 



##############
cpm_df <- read_tsv("./data/edger_cpm_fit.tsv")
edger_results <- read_tsv("./data/edger_de_genes.tsv")
#genes <- c("WDFY4", "BLK", "IKZF2", "FOSL2", "PXK", "SLC15A4", "IRF5")
genes <- c("BLK", "PXK", "WDFY4", "IKZF2")

select_genes <- edger_results %>%
    filter(gene_name %in% genes) %>%
    select(stim, gene_id, gene_name)

cpm_plot_df <- inner_join(select_genes, cpm_df) %>%
    mutate(hours = parse_number(time),
           hours = factor(hours, levels = sort(unique(hours))),
           condition = paste(stim, time, sep = "_"),
	   stim = factor(stim, levels = c("CD40L", "TLR9", "BCR", "TLR7", "BCR-TLR7", "DN2"))) %>%
    mutate_at(vars(cpm:fit), ~2^.)

p_vals <- edger_results %>%
    filter(gene_name %in% select_genes$gene_name) %>%
    select(gene_id, gene_name, stim, p = PValue, fdr = FDR) %>%
    arrange(gene_name, stim) %>%
    left_join(group_by(cpm_plot_df, gene_id, gene_name, stim) %>% summarise(cpm = max(cpm)) %>% ungroup()) %>%
    mutate(p_lab = format(p, format = "e", digits = 2),
	   p_lab = paste("p =", p_lab),
	   p_lab = ifelse(fdr < 0.05, paste(p_lab, "*"), p_lab))

plot_timecourse <- function(gene) {

    out <- cpm_plot_df %>% 
	filter(gene_name == gene) %>%
	ggplot() +
	geom_quasirandom(aes(x = hours, y = cpm, fill = condition),
			 method = "smiley", width = .2, 
			 shape = 21, stroke = .2, size = 3) +
	geom_boxplot(aes(x = hours, y = cpm, fill = condition), 
		     outlier.color = NA, size = .25, width = .5, alpha = .5) +
	scale_fill_manual(values = stim_colors) + 
	geom_text(data = p_vals %>% filter(gene_name == gene), 
		  aes(x = 0.5, y = cpm * 1.25, label = p_lab),
		   hjust = "inward", vjust = "inward", size = 4) +
	facet_grid(stim~gene_name, scale = "free_y") +
	theme(panel.grid.minor = element_blank(),
	      panel.grid.major.x = element_blank(),
	      panel.background = element_rect(fill = "grey96"),
	      legend.position = "none") +
	labs(x = NULL, y = "CPM")

    ggsave(sprintf("./plots/timecourse_%s.png", gene), out, width = 4, height = 8)
}

walk(genes, plot_timecourse)


dn2_df <- cpm_plot_df %>%
    filter(stim == "DN2") %>%
    mutate(stim = as.character(stim))

dn2_pvals <- filter(p_vals, stim == "DN2")


out <- ggplot(data = dn2_df) +
    geom_quasirandom(aes(x = hours, y = cpm, fill = condition),
		     method = "smiley", width = .2,
		     shape = 21, stroke = .2, size = 2) +
    geom_boxplot(aes(x = hours, y = cpm, fill = condition),
		 outlier.color = NA, size = .25, width = .5, alpha = .5) +
    scale_fill_manual(values = stim_colors) +
    geom_text(data = dn2_pvals, aes(x = .5, y = cpm * 1.25, label = p_lab),
	      hjust = "inward", vjust = "inward", size = 3, family = "Arial") +
    facet_wrap(~gene_name, ncol = 2, scale = "free_y") +
    theme_bw() +
    theme(text = element_text(size = 12, family = "Arial"),
	  legend.position = "none",
	  panel.grid.minor = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.grid.major.y = element_line(color = "grey96"),
	  panel.background = element_rect(fill = "white", color = "white"),
	  strip.background = element_rect(fill = "white"),
	  strip.text = element_text(face = "bold")) +
    labs(y = "CPM")

ggsave("./plots/grant.png", out, width = 3.25, height = 3.5)


