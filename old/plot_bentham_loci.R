library(tidyverse)
library(ggbeeswarm)
library(RColorBrewer)
library(readxl)
library(UpSetR)

coloc_genes <- read_lines("./colocalization/coloc_genes_bentham.txt")
coloc_genes <- coloc_genes[coloc_genes != "MHC class IIId"]
coloc_genes <- c(coloc_genes, "C4A", "C4B")

### low input
stims_order <-
    c("Unstim_0hrs", "Unstim_4hrs", "Unstim_24hrs", 
      "IL4_4hrs", "IL4_24hrs",
      paste("CD40L", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("BCR", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("TLR7", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("BCR-TLR7", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("TLR9", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("DN2", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"))
      
stim_colors <- c(brewer.pal(n = 6, "Greys")[-1],
		 brewer.pal(n = 9, "YlOrRd")[c(1, 3, 4, 5)],
		 brewer.pal(n = 9, "Blues")[c(2, 4, 6, 8)],
		 brewer.pal(n = 9, "Greens")[c(2, 4, 6, 8)],
		 grep("cyan", colors(), value = TRUE)[c(8, 2, 4, 6)],
		 grep("pink", colors(), value = TRUE)[c(16, 7, 1, 4)],
		 paste0("tomato", c("", 2:4)))

names(stim_colors) <- stims_order

sample_decode <- read_tsv("./bcell_lowinput/data/sample_decode.tsv") |>
    separate(sample_name, c("name", "stim", "time"), sep = "_", remove = FALSE) |>
    mutate(time = factor(time, levels = str_sort(unique(time), numeric = TRUE))) |>
    arrange(stim, time, name)

low_stims <- c("Unstim", "IL4", "CD40L", "TLR9", "BCR", "TLR7", "BCR-TLR7", "DN2")

lowin <- read_rds("./bcell_lowinput/data/expression.rds") |>
    filter(gene_name %in% coloc_genes) |>
    left_join(sample_decode, by = c("id" = "sample_id")) |>
    select(sample_id = name, stim, time, gene_id, gene_name, tpm) |>
    mutate(stim = factor(stim, levels = low_stims),
	   condition = paste(stim, time, sep = "_"))

lowin_tx <- read_rds("./bcell_lowinput/data/expression_transcripts.rds") |>
    filter(gene_name %in% coloc_genes) |>
    left_join(sample_decode, by = c("id" = "sample_id")) |>
    select(sample_id = name, stim, time, gene_id, gene_name, tx_id, tpm) |>
    mutate(stim = factor(stim, levels = low_stims),
	   condition = paste(stim, time, sep = "_"))

edger_gene <- read_tsv("./bcell_lowinput/data/edger_de_genes.tsv")
edger_tx <- read_tsv("./bcell_lowinput/data/edger_de_transcripts.tsv")

plot_low <- function(dat, gene) {
    filter(dat, gene_name == gene) |>
    ggplot(aes(x = time, y = tpm, fill = condition)) +
    geom_quasirandom(method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 3) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 2) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  legend.position = "none") +
    labs(title = gene)
}


# Region 1

gene_1 <- "PHTF1"

edger_gene |> filter(gene_name == gene_1, FDR < 0.05)

phtf1 <- plot_low(lowin, gene_1)

ggsave(sprintf("./plots_labmeeting/%s.png", gene_1),
       phtf1, width = 5, height = 6)




# Region 2
gene_2 <- "FCGR2A" 

fcgr2a <- plot_low(lowin, gene_2)

edger_tx |> filter(gene_name == gene_2)

gene_2.2 <- "SDHC"
plot_low(lowin, gene_2.2)

gene_2.2_tx <- edger_tx |> 
    filter(gene_name == gene_2.2, FDR < .05) |> 
    distinct(tx_id) |>
    pull(tx_id)

gene_2.2_plot <- 
    lowin_tx |> mutate(tx_id = sub("\\.\\d+$", "", tx_id)) |> 
    filter(tx_id %in% gene_2.2_tx) |>
    ggplot(aes(x = time, y = tpm, fill = condition)) +
    geom_quasirandom(method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 3) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 2) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  legend.position = "none") +
    labs(title = gene_2.2_tx)
    
ggsave(sprintf("./plots_labmeeting/%s.png", gene_2.2),
       gene_2.2_plot, width = 5, height = 6)



# Region 3
gene_3 <- "TNFSF4"
edger_gene |> filter(gene_name == gene_3, FDR < 0.05)
edger_tx |> filter(gene_name == gene_3, FDR < 0.05)


plot_3 <- plot_low(lowin, gene_3)

ggsave(sprintf("./plots_labmeeting/%s.png", gene_3),
       plot_3, width = 5, height = 6)

gene_3.2 <- "AL645568.1"
edger_gene |> filter(gene_name == gene_3.2, FDR < 0.05)
edger_tx |> filter(gene_name == gene_3.2, FDR < 0.05)


# Region 4
gene_4 <- "NCF2"
edger_gene |> filter(gene_name == gene_4, FDR < 0.05)
edger_tx |> filter(gene_name == gene_4, FDR < 0.05)

gene_4.2 <- "SMG7"
edger_gene |> filter(gene_name == gene_4.2, FDR < 0.05)
edger_tx |> filter(gene_name == gene_4.2, FDR < 0.05)


plot_4 <- plot_low(lowin, gene_4)

ggsave(sprintf("./plots_labmeeting/%s.png", gene_4),
       plot_4, width = 5, height = 6)

plot_4.2 <- plot_low(lowin, gene_4.2)

ggsave(sprintf("./plots_labmeeting/%s.png", gene_4.2),
       plot_4.2, width = 5, height = 6)



# Region 5
gene_5 <- "SPRED2"
edger_gene |> filter(gene_name == gene_5, FDR < 0.05)
edger_tx |> filter(gene_name == gene_5, FDR < 0.05)


plot_5 <- plot_low(lowin, gene_5)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_5),
       plot_5, width = 5, height = 6)


# Region 6
gene_6 <- "IFIH1"
edger_gene |> filter(gene_name == gene_6, FDR < 0.05)
edger_tx |> filter(gene_name == gene_6, FDR < 0.05)

plot_6 <- plot_low(lowin, gene_6)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_6),
       plot_6, width = 5, height = 6)

gene_6.2 <- "DPP4"
edger_gene |> filter(gene_name == gene_6.2, FDR < 0.05)
edger_tx |> filter(gene_name == gene_6.2, FDR < 0.05)

plot_6.2 <- plot_low(lowin, gene_6.2)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_6.2),
       plot_6.2, width = 5, height = 6)



# Region 7
gene_7 <- "STAT4"
edger_gene |> filter(gene_name == gene_7, FDR < 0.05)

tx_7 <- edger_tx |> 
    filter(gene_name == gene_7, FDR < 0.05) |>
    filter(F == max(F)) |> 
    pull(tx_id)

plot_7 <- plot_low(lowin, gene_7)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_7),
       plot_7, width = 5, height = 6)

plot_7_tx <- lowin_tx |> 
    mutate(tx_id = sub("\\.\\d+$", "", tx_id)) |> 
    filter(tx_id == tx_7) |>
    ggplot(aes(x = time, y = tpm, fill = condition)) +
    geom_quasirandom(method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 3) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 2) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  legend.position = "none") +
    labs(title = tx_7)

ggsave(sprintf("./plots_labmeeting/%s_tx.png", gene_7),
       plot_7_tx, width = 5, height = 6)

gene_7.2 <- "STAT1"
edger_gene |> filter(gene_name == gene_7.2, FDR < 0.05)
edger_tx |> filter(gene_name == gene_7.2, FDR < 0.05)

plot_7.2 <- plot_low(lowin, gene_7.2)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_7.2),
       plot_7.2, width = 5, height = 6)




# Region 8
gene_8 <- "IKZF2"
edger_gene |> filter(gene_name == gene_8, FDR < 0.05)
edger_tx |> filter(gene_name == gene_8, FDR < 0.05)

plot_8 <- plot_low(lowin, gene_8)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_8),
       plot_8, width = 5, height = 6)



# Region 9
gene_9 <- "PXK"
edger_gene |> filter(gene_name == gene_9, FDR < 0.05)
edger_tx |> filter(gene_name == gene_9, FDR < 0.05)

plot_9 <- plot_low(lowin, gene_9)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_9),
       plot_9, width = 5, height = 6)




# Region 10
gene_10 <- "IL12A"
edger_gene |> filter(gene_name == gene_10, FDR < 0.05)
edger_tx |> filter(gene_name == gene_10, FDR < 0.05)

plot_10 <- plot_low(lowin, gene_10)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_10),
       plot_10, width = 5, height = 6)



# Region 11
gene_11 <- "BANK1"
edger_gene |> filter(gene_name == gene_11, FDR < 0.05)
edger_tx |> filter(gene_name == gene_11, FDR < 0.05)

plot_11 <- plot_low(lowin, gene_11)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_11),
       plot_11, width = 5, height = 6)




# Region 12
gene_12 <- "TCF7"
edger_gene |> filter(gene_name == gene_12, FDR < 0.05)

edger_tx |> 
    filter(gene_name == gene_12, FDR < 0.05) |>
    filter(F == max(F))

plot_12 <- plot_low(lowin, gene_12)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_12),
       plot_12, width = 5, height = 6)


gene_12.2 <- "SKP1"
edger_gene |> filter(gene_name == gene_12.2, FDR < 0.05)
edger_tx |> filter(gene_name == gene_12.2, FDR < 0.05)

tx_12.2 <- edger_tx |> 
    filter(gene_name == gene_12.2, FDR < 0.05) |>
    filter(FDR == min(FDR)) |>
    pull(tx_id)

plot_12.2 <- plot_low(lowin, gene_12.2)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_12.2),
       plot_12.2, width = 5, height = 6)



plot_12.2_tx <- lowin_tx |> 
    mutate(tx_id = sub("\\.\\d+$", "", tx_id)) |> 
    filter(tx_id == tx_12.2) |>
    ggplot(aes(x = time, y = tpm, fill = condition)) +
    geom_quasirandom(method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 3) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 2) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  legend.position = "none") +
    labs(title = tx_12.2)

ggsave(sprintf("./plots_labmeeting/%s_tx.png", gene_12.2),
       plot_12.2_tx, width = 5, height = 6)




# Region 13
gene_13 <- "TNIP1"
edger_gene |> filter(gene_name == gene_13, FDR < 0.05)
edger_tx |> filter(gene_name == gene_13, FDR < 0.05)

plot_13 <- plot_low(lowin, gene_13)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_13),
       plot_13, width = 5, height = 6)


gene_13.2 <- "GPX3"
edger_gene |> filter(gene_name == gene_13.2, FDR < 0.05)
edger_tx |> filter(gene_name == gene_13.2, FDR < 0.05)



# Region 14
gene_14 <- "MIR3142HG"
edger_gene |> filter(gene_name == gene_14, FDR < 0.05)
edger_tx |> filter(gene_name == gene_14, FDR < 0.05)

plot_14 <- plot_low(lowin, gene_14)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_14),
       plot_14, width = 5, height = 6)




# Region 15
gene_15 <- "C4A"
edger_gene |> filter(gene_name == gene_15, FDR < 0.05)
edger_tx |> filter(gene_name == gene_15, FDR < 0.05)





# Region 16
gene_16 <- "PRDM1"
edger_gene |> filter(gene_name == gene_16, FDR < 0.05)
edger_tx |> filter(gene_name == gene_16, FDR < 0.05)

plot_16 <- plot_low(lowin, gene_16)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_16),
       plot_16, width = 5, height = 6)

gene_16.2 <- "ATG5"
edger_gene |> filter(gene_name == gene_16.2, FDR < 0.05)
edger_tx |> filter(gene_name == gene_16.2, FDR < 0.05)

plot_16.2 <- plot_low(lowin, gene_16.2)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_16.2),
       plot_16.2, width = 5, height = 6)


gene_16.3 <- "AL022067.1"
edger_gene |> filter(gene_name == gene_16.3, FDR < 0.05)
edger_tx |> filter(gene_name == gene_16.3, FDR < 0.05)




# Region 17
gene_17 <- "TNFAIP3"
edger_gene |> filter(gene_name == gene_17, FDR < 0.05)
edger_tx |> filter(gene_name == gene_17, FDR < 0.05)

plot_17 <- plot_low(lowin, gene_17)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_17),
       plot_17, width = 5, height = 6)





# Region 18
gene_18 <- "IRF5"
edger_gene |> filter(gene_name == gene_18, FDR < 0.05)
edger_tx |> 
    filter(gene_name == gene_18, FDR < 0.05) |> distinct(tx_id)
    filter(FDR == min(FDR))

plot_18 <- plot_low(lowin, gene_18)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_18),
       plot_18, width = 5, height = 6)



# Region 19
gene_19 <- "BLK"
edger_gene |> filter(gene_name == gene_19, FDR < 0.05)
edger_tx |> filter(gene_name == gene_19, FDR < 0.05)

plot_19 <- plot_low(lowin, gene_19)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_19),
       plot_19, width = 5, height = 6)


gene_19.2 <- "FAM167A"
edger_gene |> filter(gene_name == gene_19.2, FDR < 0.05)
edger_tx |> filter(gene_name == gene_19.2, FDR < 0.05)

plot_19.2 <- plot_low(lowin, gene_19.2)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_19.2),
       plot_19.2, width = 5, height = 6)



# Region 20
gene_20 <- "WDFY4"
edger_gene |> filter(gene_name == gene_20, FDR < 0.05)
edger_tx |> filter(gene_name == gene_20, FDR < 0.05)

plot_20 <- plot_low(lowin, gene_20)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_20),
       plot_20, width = 5, height = 6)

gene_20.2 <- "AC035139.1"
edger_gene |> filter(gene_name == gene_20.2, FDR < 0.05)
edger_tx |> filter(gene_name == gene_20.2, FDR < 0.05)



# Region 21
gene_21 <- "ARID5B"
edger_gene |> filter(gene_name == gene_21, FDR < 0.05)
edger_tx |> filter(gene_name == gene_21, FDR < 0.05)

plot_21 <- plot_low(lowin, gene_21)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_21),
       plot_21, width = 5, height = 6)





# Region 22
gene_22 <- "IRF7"
edger_gene |> filter(gene_name == gene_22, FDR < 0.05)
edger_tx |> filter(gene_name == gene_22, FDR < 0.05)

plot_22 <- plot_low(lowin, gene_22)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_22),
       plot_22, width = 5, height = 6)



# Region 23
gene_23 <- "CD44"
edger_gene |> filter(gene_name == gene_23, FDR < 0.05)

plot_23 <- plot_low(lowin, gene_23)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_23),
       plot_23, width = 5, height = 6)

tx_23 <- edger_tx |> 
    filter(gene_name == gene_23, FDR < 0.05) |> 
    filter(FDR == min(FDR)) |>
    pull(tx_id)

plot_23_tx <- lowin_tx |> 
    mutate(tx_id = sub("\\.\\d+$", "", tx_id)) |> 
    filter(tx_id == tx_23) |>
    ggplot(aes(x = time, y = tpm, fill = condition)) +
    geom_quasirandom(method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 3) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 2) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  legend.position = "none") +
    labs(title = tx_23)

ggsave(sprintf("./plots_labmeeting/%s_tx.png", gene_23),
       plot_23_tx, width = 5, height = 6)



# Region 24
gene_24 <- "ETS1"
edger_gene |> filter(gene_name == gene_24, FDR < 0.05)

plot_24 <- plot_low(lowin, gene_24)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_24),
       plot_24, width = 5, height = 6)

tx_24 <- edger_tx |> 
    filter(gene_name == gene_24, FDR < 0.05) |> 
    filter(FDR == min(FDR)) |>
    pull(tx_id)

plot_24_tx <- lowin_tx |> 
    mutate(tx_id = sub("\\.\\d+$", "", tx_id)) |> 
    filter(tx_id == tx_24) |>
    ggplot(aes(x = time, y = tpm, fill = condition)) +
    geom_quasirandom(method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 3) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 2) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  legend.position = "none") +
    labs(title = tx_24)

ggsave(sprintf("./plots_labmeeting/%s_tx.png", gene_24),
       plot_24_tx, width = 5, height = 6)



# Region 25
gene_25 <- "SH2B3"
edger_gene |> filter(gene_name == gene_25, FDR < 0.05)

tx_25 <- edger_tx |> 
    filter(gene_name == gene_25, FDR < 0.05) |> 
    filter(FDR == min(FDR)) |>
    pull(tx_id)

plot_25 <- plot_low(lowin, gene_25)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_25),
       plot_25, width = 5, height = 6)

plot_25_tx <- lowin_tx |> 
    mutate(tx_id = sub("\\.\\d+$", "", tx_id)) |> 
    filter(tx_id == tx_25) |>
    ggplot(aes(x = time, y = tpm, fill = condition)) +
    geom_quasirandom(method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 3) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 2) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  legend.position = "none") +
    labs(title = tx_25)

ggsave(sprintf("./plots_labmeeting/%s_tx.png", gene_25),
       plot_25_tx, width = 5, height = 6)



gene_25.2 <- "ATXN2"
edger_gene |> filter(gene_name == gene_25.2, FDR < 0.05)

tx_25.2 <- edger_tx |> 
    filter(gene_name == gene_25.2, FDR < 0.05) |> 
    filter(FDR == min(FDR)) |>
    pull(tx_id)

plot_25.2_tx <- lowin_tx |> 
    mutate(tx_id = sub("\\.\\d+$", "", tx_id)) |> 
    filter(tx_id == tx_25.2) |>
    ggplot(aes(x = time, y = tpm, fill = condition)) +
    geom_quasirandom(method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 3) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 2) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  legend.position = "none") +
    labs(title = tx_25.2)

ggsave(sprintf("./plots_labmeeting/%s_tx.png", gene_25.2),
       plot_25.2_tx, width = 5, height = 6)







# Region 26
gene_26 <- "SLC15A4"
edger_gene |> filter(gene_name == gene_26, FDR < 0.05)

edger_tx |> filter(gene_name == gene_26, FDR < 0.05) 

plot_26 <- plot_low(lowin, gene_26)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_26),
       plot_26, width = 5, height = 6)





# Region 27
gene_27 <- "CSK"
edger_gene |> filter(gene_name == gene_27, FDR < 0.05)
edger_tx |> filter(gene_name == gene_27, FDR < 0.05) 

plot_27 <- plot_low(lowin, gene_27)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_27),
       plot_27, width = 5, height = 6)





# Region 28
gene_28 <- "CLEC16A"
edger_gene |> filter(gene_name == gene_28, FDR < 0.05)
edger_tx |> filter(gene_name == gene_28, FDR < 0.05) 

gene_28.2 <- "SOCS1"
edger_gene |> filter(gene_name == gene_28.2, FDR < 0.05)
edger_tx |> filter(gene_name == gene_28.2, FDR < 0.05) 

gene_28.3 <- "CIITA"
edger_gene |> filter(gene_name == gene_28.3, FDR < 0.05)
edger_tx |> filter(gene_name == gene_28.3, FDR < 0.05) 


plot_28 <- plot_low(lowin, gene_28)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_28),
       plot_28, width = 5, height = 6)

plot_28.2 <- plot_low(lowin, gene_28.2)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_28.2),
       plot_28.2, width = 5, height = 6)





# Region 29
gene_29 <- "ITGAM"
edger_gene |> filter(gene_name == gene_29, FDR < 0.05)
edger_tx |> filter(gene_name == gene_29, FDR < 0.05) 

plot_29 <- plot_low(lowin, gene_29)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_29),
       plot_29, width = 5, height = 6)

gene_29.2 <- "ITGAX"
edger_gene |> filter(gene_name == gene_29.2, FDR < 0.05)
edger_tx |> filter(gene_name == gene_29.2, FDR < 0.05) 

plot_29.2 <- plot_low(lowin, gene_29.2)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_29.2),
       plot_29.2, width = 5, height = 6)






# Region 30
gene_30 <- "IRF8"
edger_gene |> filter(gene_name == gene_30, FDR < 0.05)
edger_tx |> filter(gene_name == gene_30, FDR < 0.05) 

plot_30 <- plot_low(lowin, gene_30)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_30),
       plot_30, width = 5, height = 6)





# Region 31
gene_31 <- "GSDMB"
edger_gene |> filter(gene_name == gene_31, FDR < 0.05)
edger_tx |> filter(gene_name == gene_31, FDR < 0.05) 

gene_31.2 <- "ORMDL3"
edger_gene |> filter(gene_name == gene_31.2, FDR < 0.05)
edger_tx |> filter(gene_name == gene_31.2, FDR < 0.05) 

gene_31.3 <- "IKZF3"
edger_gene |> filter(gene_name == gene_31.3, FDR < 0.05)
edger_tx |> filter(gene_name == gene_31.3, FDR < 0.05) 

plot_31 <- plot_low(lowin, gene_31)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_31),
       plot_31, width = 5, height = 6)

plot_31.3 <- plot_low(lowin, gene_31.3)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_31.3),
       plot_31.3, width = 5, height = 6)





# Region 32
gene_32 <- "TYK2"
edger_gene |> filter(gene_name == gene_32, FDR < 0.05)
edger_tx |> filter(gene_name == gene_32, FDR < 0.05) 

plot_32 <- plot_low(lowin, gene_32)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_32),
       plot_32, width = 5, height = 6)




# Region 33
gene_33 <- "UBE2L3"
edger_gene |> filter(gene_name == gene_33, FDR < 0.05)

tx_33 <- edger_tx |> 
    filter(gene_name == gene_33, FDR < 0.05) |>
    filter(FDR == min(FDR)) |>
    pull(tx_id)

plot_33 <- plot_low(lowin, gene_33)
ggsave(sprintf("./plots_labmeeting/%s.png", gene_33),
       plot_33, width = 5, height = 6)

plot_33_tx <- lowin_tx |> 
    mutate(tx_id = sub("\\.\\d+$", "", tx_id)) |> 
    filter(tx_id == tx_33) |>
    ggplot(aes(x = time, y = tpm, fill = condition)) +
    geom_quasirandom(method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 3) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 2) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  legend.position = "none") +
    labs(title = tx_33)

ggsave(sprintf("./plots_labmeeting/%s_tx.png", gene_33),
       plot_33_tx, width = 5, height = 6)



###


edge_list <- read_excel("/home/ch229163/Lupus.xlsx") |>
    pivot_longer(-(1:2)) |>
    filter(value == 1) |>
    select(name, subject_id = 2) |>
    {function(x) split(x, x$name)}() |>
    map("subject_id")

png("./plots_labmeeting/upset.png", res = 300, units = "in", width = 6, height = 6)
upset(fromList(edge_list), order.by = "freq", text.scale = 2)
dev.off()
