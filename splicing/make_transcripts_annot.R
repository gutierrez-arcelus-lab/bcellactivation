library(tidyverse)

annot <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf") %>%
    read_tsv(comment = "#", col_types = "c-cii-c-c",
	     col_names = c("chr", "feature", "start", "end", "strand", "info")) %>%
    mutate(chr = factor(chr, levels = str_sort(unique(chr), numeric = TRUE)))

exon_annot <- annot %>%
    filter(feature == "exon") %>%
    select(-feature) %>%
    mutate(info = str_split(info, "; "),
	   info = map(info, ~keep(., grepl("gene_id|gene_name|transcript_id", .)))) %>%
    rownames_to_column("i") %>%
    unnest(cols = c(info)) %>%
    separate(info, c("feature", "id"), sep = " ") %>%
    mutate(id = str_remove_all(id, "[\"]")) %>%
    pivot_wider(names_from = feature, values_from = id)


intron_annot <- exon_annot %>%
    group_by(gene_id) %>%
    mutate(start_gene = min(start), end_gene = max(end)) %>%
    ungroup() %>%
    arrange(chr, start_gene, start) %>%
    group_by(gene_id, gene_name, transcript_id, start_gene, end_gene) %>%
    summarise(data = bind_cols(tibble(start = c(unique(start_gene), end)),
			       tibble(end = c(start, unique(end_gene))))) %>%
    ungroup() %>%
    unnest(cols = c(data)) %>%
    mutate(start_i = ifelse(start == end | start == start_gene, start, start + 1L),
	   end_i = ifelse(start == end | end == end_gene, end, end - 1L)) %>%
    filter(start != end) %>%
    select(gene_id, gene_name, transcript_id, start = start_i, end = end_i)

annot_out <- 
    bind_rows("exon" = select(exon_annot, gene_id, gene_name, transcript_id, start, end),
              "intron" = intron_annot,
              .id = "feature") %>%
    arrange(gene_name, gene_id, transcript_id, start)

write_rds(annot_out, "./plot_data/transcripts_annot.rds")