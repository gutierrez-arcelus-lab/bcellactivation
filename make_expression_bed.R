library(tidyverse)

annotations <- 
    "./data/gencode.v38.primary_assembly.annotation.gtf.gz" %>% 
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c")

# make gene annotations as well to join with gene_bed at the end


transcript_annot <- annotations %>%
    filter(X2 == "transcript") %>%
    transmute(chr = X1,
	      gene_id = str_extract(X6, "(?<=gene_id\\s\")[^\"]+"),
	      transcript_id = str_extract(X6, "(?<=transcript_id\\s\")[^\"]+"),
	      transcript_biotype = str_extract(X6, "(?<=transcript_type\\s\")[^\"]+"),
	      start = X3, end = X4, strand = X5)

sample_info <- 
    tibble(sample_id = list.files("./bcell_quant")) %>%
    mutate(condition_id = gsub("^[^_]+_([^_]+)_([^_]+).*$", "\\1_\\2", sample_id))

quant_df <- 
    file.path("./bcell_quant", sample_info$sample_id, "quant.sf") %>% 
    setNames(sample_info$condition_id) %>%
    map_df(read_tsv, .id = "condition_id") %>%
    left_join(sample_info, by = "condition_id") %>%
    select(condition_id, transcript_id = Name, tpm = TPM)

quant_annot_df <- quant_df %>%
    left_join(transcript_annot) %>%
    select(condition_id, transcript_id, transcript_biotype, gene_id, tpm)

quant_annot_df %>%
    filter(is.na(gene_id))

gene_df <- quant_annot_df %>%
    group_by(condition_id, gene_id) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup()

expressed_genes_df <- gene_df %>%
    group_by(gene_id) %>%
    filter(any(tpm > 0)) %>%
    ungroup()

#expressed_genes <- gene_df %>%
#    group_by(gene_id) %>%
#    filter(mean(tpm>1) >= 0.5) %>%
#    ungroup()
#
gene_bed <- gene_df %>%
    pivot_wider(names_from = condition_id, values_from = tpm) %>%
    left_join(transcript_annot, by = "gene_id") %>%
    mutate(gid = gene_id) %>%
    select(`#chr` = chr, start, end, id = gene_id, gid, strd = strand, 
	   contains("hr_")) %>%
    arrange(`#chr`, start)

write_tsv(quant_annot_df, "./transcript_quants.tsv")
write_tsv(gene_bed, "./phenotypes.bed")
system("bgzip phenotypes.bed && tabix -p bed phenotypes.bed.gz")
