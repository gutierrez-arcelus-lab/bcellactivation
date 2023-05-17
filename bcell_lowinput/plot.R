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

cpm_df <- "./data/edger_cpm_fit.tsv" |>
    read_tsv() |>
    mutate(hours = parse_number(time),
           hours = factor(hours, levels = sort(unique(hours))),
           condition = paste(stim, time, sep = "_"),
	   stim = factor(stim, levels = c("CD40L", "TLR9", "BCR", "TLR7", "BCR-TLR7", "DN2"))) |>
    mutate_at(vars(cpm:fit), ~2^.)

early_genes <- cpm_df |>
    summarise(cpm = mean(cpm), .by = c(stim, time, gene_id, gene_name)) |>
    pivot_wider(names_from = time, values_from = cpm) |>
    filter(`4hrs` > `0hrs` & `4hrs` > `72hrs`) |>
    select(stim, gene_id, gene_name)

sle_genes <- read_tsv("../bcell_scrna/reported_genes.tsv")

edger_results <- read_tsv("./data/edger_de_genes.tsv")

top_sle_genes <- edger_results |>
    filter(gene_name %in% sle_genes$gene) |>
    group_by(stim) |>
    top_n(6, -log10(PValue)) |>
    ungroup() |>
    select(stim, gene_id, gene_name)

top_early_genes <- edger_results |>
    inner_join(early_genes) |>
    group_by(stim) |>
    top_n(6, -log10(PValue)) |>
    ungroup() |>
    select(stim, gene_id, gene_name)

specific_genes <- edger_results |>	
    filter(gene_name %in% sle_genes$gene,
	   logCPM > 4) |>
    select(stim, gene_id, FDR) |>
    pivot_wider(names_from = stim, values_from = FDR) |>
    filter(! CD40L < 0.1 & ! TLR9 < 0.1) |>
    pull(gene_id)

top_specific_genes <- edger_results |>
    filter(gene_id %in% specific_genes, ! stim %in% c("TLR9", "CD40L"), FDR < 0.05) |>
    group_by(stim) |>
    top_n(6, -log10(PValue)) |>
    ungroup() |>
    select(stim, gene_id, gene_name)


cpm_plot_df <- cpm_df |>
    #inner_join(top_specific_genes)
    inner_join(top_genes)

p_vals <- edger_results |>
    inner_join(distinct(cpm_plot_df, stim, gene_id, gene_name)) |>
    select(gene_id, gene_name, stim, p = PValue, fdr = FDR) |>
    arrange(gene_name, stim) |>
    left_join(summarise(cpm_plot_df, cpm = max(cpm), .by = c(gene_id, gene_name, stim))) |>
    mutate(p_lab = format(p, format = "e", digits = 2),
	   p_lab = paste("p =", p_lab),
	   p_lab = ifelse(fdr < 0.05, paste(p_lab, "*"), p_lab))



plot_tc <- function(stimulus) {
    
    cpm_plot_df |> 
	filter(stim == first(stimulus)) |>
	ggplot() +
	geom_quasirandom(aes(x = hours, y = cpm, fill = condition),
			 method = "smiley", width = .2, 
			 shape = 21, stroke = .2, size = 3) +
	geom_text(data = p_vals |> filter(stim == stimulus), 
		  aes(x = 0.5, y = cpm * 1.25, label = p_lab),
		  hjust = "inward", vjust = "inward", size = 4) +
	scale_fill_manual(values = stim_colors) + 
	facet_wrap(~gene_name, scale = "free_y", nrow = 1) +
	theme(panel.grid.minor = element_blank(),
	      panel.grid.major.x = element_blank(),
	      panel.background = element_rect(fill = "grey96"),
	      legend.position = "none") +
	labs(x = NULL, y = "CPM")
}

out_plot_list <- map(levels(cpm_df$stim), plot_tc) |>
    plot_grid(plotlist = _, ncol = 1)

out_plot_list <- map(c("BCR", "TLR7", "BCR-TLR7", "DN2"), plot_tc) |>
    plot_grid(plotlist = _, ncol = 1)

ggsave("./plots/timecourse_top_specific.png", out_plot_list, width = 10, height = 7)
ggsave("./plots/timecourse_topgenes.png", out_plot_list, width = 12, height = 10)




# IKZF1
stims_order2 <- c("Unstim_0hrs", "Unstim_4hrs", "Unstim_24hrs",
		  "IL4_4hrs", "IL4_24hrs",
      paste("CD40L", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("BCR", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("TLR7", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("BCR-TLR7", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("TLR9", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("DN2", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"))
      
stim_colors2 <- c("grey90", "grey70", "grey50",
		 "grey30", "grey10",
		 brewer.pal(n = 9, "YlOrRd")[c(1, 3, 4, 5)],
		 brewer.pal(n = 9, "Blues")[c(2, 4, 6, 8)],
		 brewer.pal(n = 9, "Greens")[c(2, 4, 6, 8)],
		 grep("cyan", colors(), value = TRUE)[c(8, 2, 4, 6)],
		 grep("pink", colors(), value = TRUE)[c(16, 7, 1, 4)],
		 paste0("tomato", c("", 2:4)))

names(stim_colors2) <- stims_order2

expr_df <- read_rds("./data/expression.rds")

meta <- read_tsv("./data/sample_decode.tsv") |>
    separate(sample_name, c("donor_id", "stim", "time"), sep = "_") |>
    unite("condition", c(stim, time), sep = "_", remove = FALSE) |>
    mutate(time = factor(time, levels = c("0hrs", "4hrs", "24hrs", "48hrs", "72hrs")),
	   stim = factor(stim, levels = c("Unstim", "IL4", "CD40L", "TLR9", "BCR", "TLR7", "BCR-TLR7", "DN2")))

expr_df_meta <- left_join(expr_df, meta, by = c("id" = "sample_id"), multiple = "all")
    


ikzf1_plot <- 
    expr_df_meta |> 
    filter(gene_name == "IKZF1") |>
    ggplot() +
	geom_quasirandom(aes(x = time, y = tpm, fill = condition),
			 method = "smiley", width = .2, 
			 shape = 21, stroke = .2, size = 3) +
	scale_fill_manual(values = stim_colors2) + 
	facet_wrap(~stim, ncol = 2) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
	      panel.grid.major.x = element_blank(),
	      panel.background = element_rect(fill = "grey96"),
	      legend.position = "none") +
	labs(x = NULL, y = "TPM")

ggsave("./plots/ikzf1.png", ikzf1_plot, width = 5, height = 5.5)

