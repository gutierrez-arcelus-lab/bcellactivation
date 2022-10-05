library(tidyverse)
library(cowplot)
library(RColorBrewer)

meta_df <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_lowinput/metadata/description_hsapiens.tsv" %>%
    read_tsv() %>%
    select(plate:time)

sample_ids <- paste(meta_df$plate, meta_df$well, sep = "_")

gene_tx <- "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.txToGene.tsv" %>%
    read_tsv()

quant_df <- file.path("./results/salmon", sample_ids, "quant.sf") %>%
    setNames(sample_ids) %>%
    map_df(read_tsv, .id = "id")

quant_gene <- quant_df %>%
    left_join(gene_tx, by = c("Name" = "transcript_id")) %>%
    group_by(id, gene_id, gene_name) %>%
    summarise(tpm = sum(TPM)) %>%
    ungroup()

variable_genes <- quant_gene %>%
    group_by(gene_id) %>%
    summarise(v = var(tpm)) %>%
    ungroup() %>%
    top_n(2000, v) %>%
    arrange(desc(v))

variable_matrix <- quant_gene %>%
    filter(gene_id %in% variable_genes$gene_id) %>%
    select(-gene_name) %>%
    pivot_wider(names_from = gene_id, values_from = tpm) %>%
    column_to_rownames("id") %>%
    as.matrix()


# PCA
pca <- prcomp(variable_matrix, center = TRUE, scale. = TRUE, rank. = 50)

pc_scores <- as_tibble(pca$x, rownames = "id")


stims_order <-
    c("Unstim_0hrs", "Unstim_4hrs", "Unstim_24hrs", 
      "IL4_4hrs", "IL4_24hrs",
      paste("CD40L", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("BCR", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("TLR7", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("BCR-TLR7", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("TLR9", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      paste("DN2", paste0(c(4, 24, 48, 72), "hrs"), sep = "_"),
      "NA_NA")
      
stim_colors <- c(brewer.pal(n = 6, "Greys")[-1],
		 paste0("goldenrod", 1:4),
		 brewer.pal(n = 5, "Blues")[-1],
		 brewer.pal(n = 5, "Greens")[-1],
		 paste0("cyan", 1:4),
		 paste0("mediumpurple", 1:4),
		 paste0("tomato", 1:4),
		 "white")

names(stim_colors) <- stims_order


pca_df <- pc_scores %>%
    select(id, PC1:PC8) %>%
    separate(id, c("plate", "well"), sep = "_") %>%
    left_join(meta_df, by = c("plate", "well")) %>%
    unite(stim, c("stim", "time"), sep = "_") %>%
    mutate(stim = factor(stim, levels = stims_order))

pc1.2 <- ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(aes(fill = stim), size = 4, shape = 21) +
    scale_fill_manual(values = stim_colors) +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  legend.key.height = unit(.5, "cm")) +
    guides(fill = guide_legend(ncol = 1))

pc3.4 <- ggplot(pca_df, aes(PC3, PC4)) +
    geom_point(aes(fill = stim), size = 4, shape = 21) +
    scale_fill_manual(values = stim_colors) +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  legend.key.height = unit(.5, "cm")) +
    guides(fill = guide_legend(ncol = 1))

pc5.6 <- ggplot(pca_df, aes(PC5, PC6)) +
    geom_point(aes(fill = stim), size = 4, shape = 21) +
    scale_fill_manual(values = stim_colors) +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  legend.key.height = unit(.5, "cm")) +
    guides(fill = guide_legend(ncol = 1))

pc7.8 <- ggplot(pca_df, aes(PC7, PC8)) +
    geom_point(aes(fill = stim), size = 4, shape = 21) +
    scale_fill_manual(values = stim_colors) +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  legend.key.height = unit(.5, "cm")) +
    guides(fill = guide_legend(ncol = 1))


plot_grid(plot_grid(pc1.2 + theme(legend.position = "none"), 
		    pc3.4 + theme(legend.position = "none"), 
		    pc5.6 + theme(legend.position = "none"), 
		    pc7.8 + theme(legend.position = "none"),
		    ncol = 2),
	  get_legend(pc1.2),
	  nrow = 1, rel_widths = c(1, .25))


# UMAP
library(uwot)
set.seed(1)
umap_df <- pc_scores %>%
    column_to_rownames("id") %>%
    select(PC1:PC10) %>%
    umap() %>%
    as_tibble(rownames = "id")

umap_df %>%
    separate(id, c("plate", "well"), sep = "_") %>%
    left_join(meta_df, by = c("plate", "well")) %>%
    unite(stim, c("stim", "time"), sep = "_") %>%
    mutate(stim = factor(stim, levels = stims_order)) %>%
    ggplot(aes(V1, V2)) +
    geom_point(aes(fill = stim), size = 4, shape = 21) +
    scale_fill_manual(values = stim_colors) +
    theme_bw()


# FASTQC plots
meta_long <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_lowinput/metadata/description_hsapiens_longformat.tsv" %>%
    read_tsv() %>%
    pivot_longer(fq1:fq2, names_to = "dummy", values_to = "fastq") %>%
    mutate(fastq = basename(fastq),
	   fastq = sub("\\.fastq\\.gz", "", fastq)) %>%
    select(-barcode_seq, -dummy)

fastqc_1 <- 
    "/home/ch229163/labshr/Lab_datasets/B_cells_lowinput/SK-56QR/get.broadinstitute.org/pkgs/SN0263576/multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt" %>%
    read_tsv()

fastqc_2 <- 
    "/home/ch229163/labshr/Lab_datasets/B_cells_lowinput/SK-56QS/get.broadinstitute.org/pkgs/SN0263542/multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt" %>%
    read_tsv()

fastqc_df <- bind_rows(fastqc_1, fastqc_2) %>%
    left_join(meta_long, ., c("fastq" = "Sample")) %>%
    pivot_longer(`Unique Reads`:`Duplicate Reads`, names_to = "read_type", values_to = "n") %>%
    mutate(read_type = sub("\\sReads$", "", read_type),
	   fastq = sub("^.+([12])$", "\\1", fastq),
	   id = paste(sample_id, stim, time, paste0("L", lane), paste0("fq", fastq), sep = "_"))

fastqc_df %>%
    group_by(sample_id, stim, time) %>%
    mutate(mn = min(n[read_type == "Unique"])) %>%
    ungroup() %>%
    arrange(plate, mn) %>%
    mutate(id = fct_inorder(id)) %>%
    ggplot(aes(id, n, fill = read_type)) +
    geom_col(width = 1.05, show.legend = FALSE) +
    ggsci::scale_fill_npg() +
    facet_wrap(~plate, ncol = 1, scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 2.5, hjust = 1, vjust = 1, angle = 90),
	  panel.grid = element_blank()) +
    labs(x = NULL, y = "Total reads")

ggsave("./seqdepth.png", width = 20)

fastqc_df %>%
    filter(fastq == 1, grepl("rep", sample_id)) %>%
    distinct(sample_id, stim, time)

fastqc_df %>%
    filter(fastq == 1, !grepl("rep", sample_id), sample_id != "BLANK") %>%
    group_by(plate, sample_id, stim, time) %>%
    summarise(n = sum(n)) %>%
    ungroup() %>%
    unite("stim", c("stim", "time"), sep = "_") %>%
    mutate(stim = factor(stim, levels = stims_order[-length(stims_order)])) %>%
    ggplot(aes(sample_id, n)) +
    geom_col(aes(fill = stim), position = "dodge", 
	     color = "black", size = .1) +
    scale_fill_manual(values = stim_colors[-length(stim_colors)]) +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(~plate, scales = "free_x", ncol = 1) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = NULL, y = "Total reads")

ggsave("seqdepth_conditions.png", width = 10, height = 5)


order2 <- stims_order %>%
    grep("NA", ., invert = TRUE, value = TRUE) %>%
    sub("^([^_]+)_.+$", "\\1", .) %>%
    unique()

fastqc_df %>%
    filter(fastq == 1, !grepl("rep", sample_id), sample_id != "BLANK") %>%
    group_by(sample_id, stim, time) %>%
    summarise(n = sum(n)) %>%
    ungroup() %>%
    unite("stim_", c("stim", "time"), sep = "_", remove = FALSE) %>%
    mutate(stim = factor(stim, levels = order2),
	   stim_ = factor(stim_, levels = stims_order[-length(stims_order)]),
	   time = factor(time, levels = str_sort(unique(time), numeric = TRUE))) %>%
    ggplot(aes(sample_id, n)) +
    geom_col(aes(fill = stim_), position = "dodge", color = "black", size = .1) +
    scale_fill_manual(values = stim_colors[-length(stim_colors)]) +
    scale_y_continuous(labels = scales::comma) +
    facet_grid(time~stim) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(n = NULL, y = "Total reads", fill = "Stim")
   
ggsave("seqdepth_stimvstime.png")

fastqc_df %>%
    filter(fastq == 1) %>%
    group_by(plate, well, sample_id, stim, time) %>%
    summarise(n = sum(n)) %>%
    ungroup() %>%
    extract(well, c("row", "column"), "([A-H])(\\d+)", convert = TRUE) %>%
    arrange(plate, row, column) %>%
    unite(id, c(sample_id, stim, time), sep = "\n") %>%
    mutate(row = factor(row, levels = LETTERS[8:1]),
	   column = factor(column)) %>%
    ggplot(data = ., aes(column, row)) +
	geom_tile(aes(fill = n), color = "black") +
	geom_text(aes(label = id), size = 2.5, lineheight = .8) +
	scale_fill_gradient(low = "white", high = "tomato4",
			     labels = scales::comma) +
	facet_wrap(~plate, ncol = 1) +
	theme_minimal() +
	theme(axis.title = element_blank(),
	      plot.margin = margin(t = 1, b = 1, r = 1, l = 1, unit = "cm"),
	      plot.background = element_rect(fill = "white", color = "white")) +
    labs(fill = "Total\nreads")

ggsave("seqdepth_plate.png"
