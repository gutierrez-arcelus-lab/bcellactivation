library(tidyverse)
library(ggsci)
library(ggnewscale)

bams <- read_lines("./bam_list.txt")

bam_df <- tibble(path = bams) %>%
    separate(path, c("dummy", "dataset", "dir", "bam"), sep = "/") %>%
    mutate(depth = sub("Aligned.sortedByCoord.out.bam", "depth.txt", bam)) %>%
    unite(path, c(dummy, dataset, dir, depth), sep = "/", remove = FALSE) %>%
    mutate(covdf = map(path, ~read_tsv(., col_names = FALSE)) )

bam_df %>% select(dataset, bam, covdf)

txi <- read_rds("../splicing/plot_data/transcripts_annot.rds")

txi_y <- txi %>% 
    filter(gene_name == "EIF4A2") %>%
    distinct(transcript_id) %>%
    mutate(y = seq(5000, by = 250, length.out = n()))

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
    unnest(covdf)

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

ggsave("./eif4a2.png", out, width = 10, height = 6)

