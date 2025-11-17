library(tidyverse)

# Gencode annotations *Need RAM
annotations <- 
    "../data/gencode.v38.primary_assembly.annotation.gtf" %>% 
    read_tsv(comment = "#", col_types = "c-cii-c-c",
	     col_names = c("chr", "feature", "start", "end", "strand", "info"))

# Gene annotations
gene_annot <- annotations %>%
    filter(feature == "gene") %>%
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+"),
	   gene_type = str_extract(info, "(?<=gene_type\\s\")[^\"]+")) %>%
    select(-info)

# Transcript annotations
transcript_annot <- annotations %>%
    filter(feature == "transcript") %>%
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   transcript_id = str_extract(info, "(?<=transcript_id\\s\")[^\"]+"),
	   transcript_type = str_extract(info, "(?<=transcript_type\\s\")[^\"]+")) %>%
    select(-info)

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

# Join annotations
quant_annot_df <- quant_df %>%
    left_join(transcript_annot) %>%
    select(condition_id, transcript_id, transcript_type, gene_id, tpm)

# Gene-level estimates
gene_df <- quant_annot_df %>%
    group_by(condition_id, gene_id) %>%
    summarise(tpm = sum(tpm)) %>%
    ungroup()

# Select expressed genes
expressed_genes_df <- gene_df %>%
    group_by(gene_id) %>%
    filter(any(tpm > 10)) %>%
    ungroup()

expressed_genes_df %>%
    distinct(gene_id)

# Log Fold change
logfc_df <- expressed_genes_df %>%
    mutate(tpm = tpm + 1L) %>%
    pivot_wider(names_from = condition_id, values_from = tpm) %>%
    pivot_longer(-(1:2), names_to = "condition_id", values_to = "tpm") %>%
    select(gene_id, resting = 2, condition_id, tpm) %>%
    mutate(log2fc = log2(tpm/resting)) %>%
    left_join(gene_annot, by = "gene_id") %>%
    select(gene_id, gene_type, condition_id, tpm, log2fc)


# BED file of gene quantifications
gene_bed <- gene_df %>%
    pivot_wider(names_from = condition_id, values_from = tpm) %>%
    left_join(gene_annot, by = "gene_id") %>%
    mutate(gid = gene_id) %>%
    select(`#chr` = chr, start, end, id = gene_id, gid, strd = strand, 
	   contains("hr_")) %>%
    arrange(`#chr`, start)

# Write output files
write_tsv(quant_annot_df, "./results/transcript_quants.tsv")
write_tsv(gene_df, "./results/gene_quants.tsv")
write_tsv(logfc_df, "./results/logfc.tsv")

write_tsv(gene_bed, "./results/phenotypes.bed")
system("bgzip ./results/phenotypes.bed && tabix -p bed ./results/phenotypes.bed.gz")
