library(tidyverse)

annotations <- 
    "./data/gencode.v38.primary_assembly.annotation.gtf.gz" %>% 
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c")

transcript_annot <- annotations %>%
    filter(X2 == "transcript") %>%
    transmute(chr = X1,
	      gene_id = str_extract(X6, "(?<=gene_id\\s\")[^\"]+"),
	      gene_name = str_extract(X6, "(?<=gene_name\\s\")[^\"]+"),
	      transcript_id = str_extract(X6, "(?<=transcript_id\\s\")[^\"]+"),
	      start = X3, end = X4, strand = X5)

sample_info <- 
    tibble(sample_id = list.files("./bcell_quant")) %>%
    mutate(id = gsub("^[^_]+_([^_]+)_([^_]+).*$", "\\1_\\2", sample_id))

quant_df <- 
    file.path("./bcell_quant", sample_info$sample_id, "quant.sf") %>% 
    setNames(sample_info$sample_id) %>%
    map_df(read_tsv, .id = "sample_id") %>%
    left_join(sample_info, by = "sample_id") %>%
    select(id, transcript_id = Name, tpm = TPM)

gene_df <- quant_df %>%
    left_join(select(transcript_annot, transcript_id, gene_id, gene_name)) %>%
    group_by(id, gene_id) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup()

expressed_genes <- gene_df %>%
    group_by(gene_id) %>%
    filter(mean(tpm>1) >= 0.5) %>%
    ungroup()

gene_bed <- expressed_genes %>%
    pivot_wider(names_from = id, values_from = tpm) %>%
    left_join(transcript_annot, by = "gene_id") %>%
    mutate(gid = gene_id) %>%
    select(`#chr` = chr, start, end, id = gene_id, gid, strd = strand, 
	   contains("hr_")) %>%
    arrange(`#chr`, start)

write_tsv(gene_bed, "phenotypes.bed")
system("bgzip phenotypes.bed && tabix -p bed phenotypes.bed.gz")
