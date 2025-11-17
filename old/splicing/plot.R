library(tidyverse)

annot <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf.gz" |>
    read_tsv(comment = "##", col_names = FALSE) |> 
    filter(X3 == "transcript") |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	   transcript_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+")) |> 
    select(gene_id, gene_name, transcript_id)    

scharer_info <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/Scharer/RNAseq/metadata/SraRunTable.txt" %>%
    read_csv() %>%
    select(run = Run, sample_name = `Sample Name`,  cell_subtype, subject_status)

expr_df <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/read_mapping/scharer/salmon", 
	      scharer_info$run, 
	      "quant.sf") %>%
    setNames(scharer_info$run) %>%
    map_df(read_tsv, .id = "sampleid")

expr_annot_df <- expr_df %>%
    inner_join(annot, by = c("Name" = "transcript_id")) %>%
    left_join(select(scharer_info, sampleid = run, status = subject_status, cell_subtype), by = "sampleid") %>%
    select(sampleid, status, cell_subtype, gene_id, gene_name, transcript_id = Name, counts = NumReads, tpm = TPM)


dn2_samples <- filter(scharer_info, grepl("DN2", cell_subtype)) %>%
    select(sampleid = run, status = subject_status)

irf5 <- filter(expr_annot_df, gene_name == "IRF5", sampleid %in% dn2_samples$sampleid)

txiannots <- read_rds("./plot_data/transcripts_annot.rds") %>%
    mutate(transcript_id = str_remove(transcript_id, "\\.\\d+$"))

load("./cluster_perdataset/results/scharer/DN/DN.Rdata")

intron_df <- as_tibble(introns) %>%
    filter(gene == "IRF5")

plot_df <- txiannots %>% 
    filter(gene_name == "IRF5") %>%
    group_by(transcript_id) %>%
    mutate(i = row_number(),
           col = case_when(feature == "intron" & i == min(i) ~ NA_character_,
                           feature == "intron" & i == max(i) ~ NA_character_,
                           feature == "exon" ~ "exon",
                           TRUE ~ "intron")) %>%
    ungroup()

irf5_struc_plot <- plot_df %>%
    filter(!is.na(col)) %>%
    filter(start > 128946500 & end < 128948000) %>%
    ggplot() +
    geom_vline(xintercept = unique(intron_df$start), linetype = 2) +
    geom_vline(xintercept = intron_df$end, linetype = 2) +
    geom_segment(aes(x = start, xend = end, 
                     y = transcript_id, yend = transcript_id,
                     color = as.factor(col)), size = 5) +
    scale_x_continuous(labels = function(x) x/1e6L,
		       breaks = scales::pretty_breaks(6)) +
    scale_color_manual(values = c("intron" = "grey80", "exon" = "blue"),
		       na.translate = FALSE) +
    facet_wrap(~gene_name, scales = "free", ncol = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
	  panel.background = element_rect(color = "white", fill = "white")) +
    labs(x = "Position (Mb)", y = NULL, color = "Feature")

ggsave("./plots/irf5_structure.png", irf5_struc_plot)


irf5_exp_plot <- irf5 %>%
    mutate(status = ifelse(grepl("SLE", status), "SLE", status)) %>%
    mutate(transcript_id = str_remove(transcript_id, "\\.\\d+$")) %>%
    ggplot(aes(x = reorder(transcript_id, tpm, "median"), 
	       y = tpm)) +
    geom_boxplot(aes(fill = status, color = status), outlier.color = NA, alpha = .25) +
    geom_jitter(aes(group = status, color = status), 
		size = 2, 
		position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c("Healthy" = "cornflowerblue", "SLE" = "tomato3")) +
    scale_fill_manual(values = c("Healthy" = "cornflowerblue", "SLE" = "tomato3")) +
    facet_wrap(~gene_name, scales = "free", ncol = 3) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
	  text = element_text(size = 12),
	  axis.text.x = element_text(angle = 90)) +
    labs(x = NULL, y = "TPM")

ggsave("./plots/irf5_expression.png", irf5_exp_plot)


cpm_df <- read_tsv("../bcell_lowinput/data/edger_cpm_fit.tsv")

cpm_irf <- filter(cpm_df, gene_name == "IRF5") %>%
    mutate(hours = parse_number(time),
           hours = factor(hours, levels = sort(unique(hours))),
           condition = paste(stim, time, sep = "_"))

library(cowplot)
library(ggbeeswarm)
library(RColorBrewer)

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
		 brewer.pal(n = 9, "YlOrRd")[c(1, 3, 4, 5)],
		 brewer.pal(n = 9, "Blues")[c(2, 4, 6, 8)],
		 brewer.pal(n = 9, "Greens")[c(2, 4, 6, 8)],
		 grep("cyan", colors(), value = TRUE)[c(8, 2, 4, 6)],
		 grep("pink", colors(), value = TRUE)[c(16, 7, 1, 4)],
		 paste0("tomato", c("", 2:4)),
		 "white")

names(stim_colors) <- stims_order

irf5_stims_plot <- cpm_irf %>%
    split(.$stim) %>%
    map(~ggplot(., aes(x = hours, y = cpm, fill = condition)) +
	    geom_jitter(shape = 21, stroke = .2, size = 3, width = .2) +
            scale_fill_manual(values = stim_colors) + 
            facet_grid(stim~gene_name) +
            theme_bw() +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major.x = element_blank(),
                  legend.position = "none") +
            labs(x = NULL)) %>%
    plot_grid(plotlist = ., ncol = 1)

ggsave("./plots/irf5_stims.png", irf5_stims_plot, width = 4, height = 8)

#IKZF2
helios_df <- 
    expr_annot_df |>
    filter(gene_name == "IKZF2") |>
    mutate(transcript_id = str_remove(transcript_id, "\\.\\d+$")) |>
    select(-gene_id, -gene_name)

helios_plot <- 
    ggplot(helios_df, aes(x = tpm, y = transcript_id)) +
    geom_boxplot(aes(fill = status, color = status), outlier.color = NA, alpha = .25) +
    geom_jitter(aes(group = status, color = status), 
		size = 1, position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c("midnightblue", "firebrick")) +
    scale_fill_manual(values = c("midnightblue", "firebrick")) +
    facet_wrap(~cell_subtype, nrow = 1, labeller = labeller(cell_subtype = label_wrap_gen(10))) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
	  legend.position = "top") +
    labs(x = "Transcripts per Million", y = NULL)

ggsave("./plots/ikzf2.png", helios_plot, height = 5, width = 8)









