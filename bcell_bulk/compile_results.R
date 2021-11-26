library(tidyverse)

# Gencode annotations
annotations <- 
    "../data/gencode.v38.primary_assembly.annotation.gtf" %>% 
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c") %>%
    mutate(X6 = trimws(X6))

# Gene annotations
gene_annot <- annotations %>%
    filter(X2 == "gene") %>%
    mutate(gene_id = str_extract(X6, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X6, "(?<=gene_name\\s\")[^\"]+"),
	   gene_type = str_extract(X6, "(?<=gene_type\\s\")[^\"]+")) %>%
    select(chr = X1, start = X3, end = X4, strand = X5, gene_id, gene_name, gene_type)


gene_annot %>%
    filter(grepl("ENS", gene_name))

# Transcript annotations
transcript_annot <- annotations %>%
    filter(X2 == "transcript") %>%
    mutate(gene_id = str_extract(X6, "(?<=gene_id\\s\")[^\"]+"),
	   transcript_id = str_extract(X6, "(?<=transcript_id\\s\")[^\"]+"),
	   transcript_type = str_extract(X6, "(?<=transcript_type\\s\")[^\"]+")) %>%
    select(chr = X1, start = X3, end = X4, strand = X5, 
	   gene_id, transcript_id, transcript_type)

select(transcript_annot, target_id = transcript_id, gene_id) %>%
    left_join(select(gene_annot, gene_id, gene_name, gene_type), by = "gene_id") %>%
    write_tsv("./data/transc_to_gene.tsv")

# Sample annotations
sample_info <- 
    tibble(sample_id = list.files("results/salmon")) %>%
    mutate(condition_id = gsub("^[^_]+_([^_]+)_([^_]+).*$", "\\1_\\2", sample_id))

# Expression estimates
quant_df <- 
    file.path("./results/salmon", sample_info$sample_id, "quant.sf") %>% 
    setNames(sample_info$condition_id) %>%
    map_df(read_tsv, .id = "condition_id") %>%
    left_join(sample_info, by = "condition_id") %>%
    select(condition_id, transcript_id = Name, counts = NumReads, tpm = TPM)

# Join annotations, compute CPM
quant_annot_df <- quant_df %>%
    group_by(condition_id) %>%
    mutate(cpm = counts/sum(counts) * 1e6) %>%
    ungroup() %>%
    left_join(transcript_annot) %>%
    select(condition_id, transcript_id, transcript_type, gene_id, cpm, tpm)

# Gene-level estimates
gene_df <- quant_annot_df %>%
    group_by(condition_id, gene_id) %>%
    summarise(cpm = sum(cpm),
	      tpm = sum(tpm)) %>%
    ungroup()

# Select expressed genes
expressed_genes_df <- gene_df %>%
    group_by(gene_id) %>%
    filter(any(tpm > 10)) %>%
    ungroup()

expressed_genes_df %>%
    distinct(gene_id)

# or:
#expressed_genes <- gene_df %>%
#    group_by(gene_id) %>%
#    filter(mean(tpm>1) >= 0.5) %>%
#    ungroup()
#

# Log Fold change
logfc_df <- expressed_genes_df %>%
    select(-tpm) %>%
    pivot_wider(names_from = condition_id, values_from = cpm) %>%
    pivot_longer(-(1:2), names_to = "condition_id", values_to = "cpm") %>%
    select(gene_id, resting = 2, condition_id, cpm) %>%
    mutate(cpm = cpm + 1L,
	   log2fc = log2(cpm/resting)) %>%
    left_join(gene_annot, by = "gene_id") %>%
    select(gene_id, gene_type, condition_id, cpm, log2fc)


# BED file of gene quantifications
gene_bed <- gene_df %>%
    select(-cpm) %>%
    pivot_wider(names_from = condition_id, values_from = tpm) %>%
    left_join(gene_annot, by = "gene_id") %>%
    mutate(gid = gene_id) %>%
    select(`#chr` = chr, start, end, id = gene_id, gid, strd = strand, 
	   contains("hr_")) %>%
    arrange(`#chr`, start)

# Write output files
write_tsv(quant_annot_df, "./data/transcript_quants.tsv")
write_tsv(gene_df, "./data/gene_quants.tsv")
write_tsv(logfc_df, "./results/logfc.tsv")

write_tsv(gene_bed, "./data/phenotypes.bed")
system("bgzip ./data/phenotypes.bed && tabix -p bed ./data/phenotypes.bed.gz")
