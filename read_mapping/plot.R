library(tidyverse)

scharer_meta <- "/lab-share/IM-Gutierrez-e2/Public/External_datasets/Scharer/RNAseq/metadata/file_description.tsv" %>%
    read_tsv() %>%
    select(sampleid = run, type = cell_subtype, status) %>%
    mutate(type = str_extract(type, "\\(.+\\)"),
	   type = gsub("[\\(\\)]", "", type))



bams <- read_lines("./bam_list.txt")

bam_df <- tibble(path = bams) %>%
    separate(path, c("dummy", "dataset", "dir", "bam"), sep = "/") %>%
    mutate(depth = sub("Aligned.sortedByCoord.out.bam", "depth.txt", bam)) %>%
    unite(path, c(dummy, dataset, dir, depth), sep = "/", remove = FALSE) %>%
    mutate(covdf = map(path, ~read_tsv(., col_names = FALSE)) )

txi <- read_rds("../splicing/plot_data/transcripts_annot.rds")

txi_y <- txi %>% 
    filter(gene_name == "EIF4A2") %>%
    distinct(transcript_id) %>%
    mutate(y = seq(5000, by = 250, length.out = n()),
	   ystd = seq(1.05, by = .05, length.out = n()))

txi_df <- txi %>% 
    filter(gene_name == "EIF4A2") %>%
    group_by(transcript_id) %>%
    mutate(i = row_number(),
           col = case_when(feature == "intron" & i == min(i) ~ NA_character_,
                           feature == "intron" & i == max(i) ~ NA_character_,
                           feature == "exon" ~ "exon",
                           TRUE ~ "intron")) %>%
    ungroup() %>%
    filter(!is.na(col)) %>%
    left_join(txi_y) %>%
    mutate(transcript_id = sub("\\.\\d+$", "", transcript_id)) %>%
    select(-i)

plot_df <- bam_df %>%
    select(dataset, bam, covdf) %>%
    mutate(sampleid = sub("^([^_]+).+$", "\\1", bam)) %>%
    select(dataset, sampleid, covdf) %>%
    unnest(covdf) %>%
    group_by(dataset, sampleid) %>%
    mutate(stdcov = X3/max(X3)) %>%
    ungroup()

out <- ggplot() +
    geom_line(data = plot_df,
	      aes(x = X2, y = X3, group = sampleid, color = dataset), 
	      alpha = .5) +
    scale_color_viridis_d() +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(breaks = seq(0, 5000, 1000)) +
    geom_segment(data = txi_df %>% 
		 filter(feature == "exon"), 
		 aes(x = start, xend = end, 
		     y = y, yend = y),
		 color = "midnightblue", size = 2.5) +
    geom_segment(data = txi_df %>% 
		 filter(feature == "intron"), 
		 aes(x = start, xend = end, 
		     y = y, yend = y),
		 color = "midnightblue", size = .5) +
    theme(panel.background = element_rect(fill = "grey98"),
	  legend.position = "bottom",
	  axis.title.y = element_text(hjust = 0.25),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank()) +
    labs(x = "Position in chr3", y = "Depth", color = "Dataset:",
	 title = expression(paste("Read coverage at the ", italic(" EIF4A2 "), " gene")),
	 subtitle = "Gencode transcripts v37") +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))


out2 <- ggplot() +
    geom_line(data = plot_df,
	      aes(x = X2, y = stdcov, group = sampleid, color = dataset), 
	      alpha = .25) +
    scale_color_viridis_d() +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    geom_segment(data = txi_df %>% 
		 filter(feature == "exon"), 
		 aes(x = start, xend = end, 
		     y = ystd, yend = ystd),
		 color = "midnightblue", size = 2.5) +
    geom_segment(data = txi_df %>% 
		 filter(feature == "intron"), 
		 aes(x = start, xend = end, 
		     y = ystd, yend = ystd),
		 color = "midnightblue", size = .5) +
    geom_text(data = distinct(txi_df, transcript_id, y, ystd),
	      aes(x = min(txi_df$start), y = ystd, label = transcript_id),
	      size = 3, color = "grey40", hjust = "outward") +
    coord_cartesian(clip = 'off') +
    theme(panel.background = element_rect(fill = "grey98"),
	  legend.position = "bottom",
	  axis.title.y = element_text(hjust = 0.25),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank()) +
    labs(x = "Position in chr3", y = "Depth", color = "Dataset:",
	 title = expression(paste("Read coverage at the ", italic(" EIF4A2 "), " gene")),
	 subtitle = "Gencode transcripts v37") +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))

ggsave("./eif4a2.png", out, width = 10, height = 6)
ggsave("./eif4a2_std.png", out2, width = 10, height = 6)


ggplot() +
    geom_line(data = plot_df %>% 
	      inner_join(scharer_meta) %>%
	      filter(type == "DN2", status == "Healthy"),
	      aes(x = X2, y = stdcov, group = sampleid, color = dataset), 
	      alpha = .25) +
    scale_color_viridis_d() +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    geom_segment(data = txi_df %>% 
		 filter(feature == "exon"), 
		 aes(x = start, xend = end, 
		     y = ystd, yend = ystd),
		 color = "midnightblue", size = 2.5) +
    geom_segment(data = txi_df %>% 
		 filter(feature == "intron"), 
		 aes(x = start, xend = end, 
		     y = ystd, yend = ystd),
		 color = "midnightblue", size = .5) +
    geom_text(data = distinct(txi_df, transcript_id, y, ystd),
	      aes(x = min(txi_df$start), y = ystd, label = transcript_id),
	      size = 3, color = "grey40", hjust = "outward") +
    coord_cartesian(clip = 'off') +
    theme(panel.background = element_rect(fill = "grey98"),
	  legend.position = "bottom",
	  axis.title.y = element_text(hjust = 0.25),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank()) +
    labs(x = "Position in chr3", y = "Depth", color = "Dataset:",
	 title = expression(paste("Read coverage at the ", italic(" EIF4A2 "), " gene")),
	 subtitle = "Gencode transcripts v37") +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))



transcript_annot <- read_tsv("../data/transcript_annots_v38.tsv")

scharer_tpm <- 
    "../read_mapping/scharer/compiled_expression_humangenes.tsv" %>%
    read_tsv() %>%
    pivot_longer(-transcript_id, names_to = "sampleid", values_to = "tpm") %>%
    inner_join(transcript_annot, by = "transcript_id") %>%
    left_join(scharer_meta, by = "sampleid") %>%
    select(sampleid, type, status, gene_id, gene_name, transcript_id, transcript_type, tpm) %>%
    mutate(transcript_id = str_remove(transcript_id, "\\.\\d+$"),
	   status = ifelse(grepl("SLE", status), "SLE", status))
    
scharer_plot_df <- scharer_tpm %>%
    filter(type == "DN2") %>%
    group_by(sampleid, type, status, transcript_type) %>%
    summarise(tpm = sum(tpm)) %>%
    group_by(sampleid, type, status) %>%
    mutate(prop = tpm/sum(tpm)) %>%
    group_by(type, status, transcript_type) %>%
    summarise(m = mean(prop)) %>%
    ungroup() %>%
    arrange(status, desc(m)) %>%
    group_by(transcript_type) %>%
    mutate(i = ifelse(any(m > 0.01), transcript_type, "other")) %>%
    ungroup() %>%
    select(type, status, transcript_type = i, m) %>%
    arrange(type, status, desc(m)) %>%
    mutate(transcript_type = fct_inorder(transcript_type),
	   transcript_type = factor(transcript_type, levels = rev(levels(transcript_type))),
	   transcript_type = fct_relevel(transcript_type, "other", after = 0))

library(ggsci)



out3 <- scharer_plot_df %>%
    ggplot(aes(x = status, y = m, fill = transcript_type)) +
    geom_col() +
    scale_fill_manual(values = c("grey", pal_npg()(9))) +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    theme(plot.background = element_rect(color = "white", fill = "white"),
	  panel.grid = element_blank()) +
    labs(y = "Average % of total expression",
	 title = "Transcriptome composition in DN2 cells\nfrom Scharer et al. dataset")

ggsave("./transcriptome_composition.png", out3)


