library(tidyverse)

# Europeans
samples_phase3 <- 
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" %>%
    read_tsv(skip = 1, col_names = FALSE) %>%
    select(subject = X1, pop = X2, continent = X3, sex = X4)

geuvadis_info <- 
    "http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt" %>%
    read_tsv() %>% 
    select(kgp_id = `Source Name`, 
           sampleid = `Comment[ENA_RUN]`,
           lab = Performer) %>% 
    distinct() %>%
    inner_join(samples_phase3, by = c("kgp_id" = "subject"))

eur_geuvadis <- filter(geuvadis_info, pop != "YRI") %>%
    select(kgp_id, sampleid)

# EBV sequence annotation
ebv_seqs <- Biostrings::readDNAStringSet("../data/ebv_genbank.fasta")

protein_names <- names(ebv_seqs) %>%
    str_remove("^lcl\\|") %>%
    str_extract("(?<=protein=).+") %>%
    str_remove("^protein ") %>%
    str_remove("\\].*$")

seq_names <- names(ebv_seqs) %>%
    str_remove("^lcl\\|") %>%
    str_remove("\\s+.*$")

gene_names <- names(ebv_seqs) %>%
    str_extract("(?<=gene=)[^\\]]+")

ebv_seq_df <- tibble(id = seq_names, gene = gene_names, protein = protein_names)

# Expression levels
quants <- read_tsv("./geuvadis_salmon_quants_ebv.bed") %>%
    select(chr = `#chr`, id, gid, starts_with("ERR")) %>%
    pivot_longer(starts_with("ERR"), names_to = "sampleid", values_to = "tpm") %>%
    inner_join(eur_geuvadis, by = "sampleid") %>%
    select(sampleid = kgp_id, chr, id, gid, tpm)

expressed_genes <- quants %>%
    group_by(id, gid) %>%
    filter(mean(tpm > 5) >= 0.75) %>%
    mutate(normexp = GenABEL::rntransform(tpm)) %>%
    ungroup()

human <- filter(expressed_genes, grepl("^chr", chr))

ebv <- filter(expressed_genes, chr == "ebv") %>%
    mutate(id = str_extract(id, "NC.+$"))

cor_df <- ebv %>%
    split(.$id) %>%
    map_df(~left_join(., human, by = "sampleid") %>%
	   group_by(id.y, gid.y) %>%
	   summarise(corr = broom::tidy(cor.test(normexp.x, normexp.y, method = "pearson", exact = FALSE))) %>%
	   ungroup() %>%
	   unnest(cols = c(corr)) %>%
	   select(gene_id_human = id.y, estimate, p = p.value), .id = "id") 

cor_annot_df <- left_join(cor_df, distinct(human, gene_id_human = id, gene_name_human = gid)) %>%
    left_join(ebv_seq_df, by = c("id")) %>%
    select(id, gene_id_human, gene_name_human, gene_name_ebv = gene, protein_name_ebv = protein, estimate, p)

#write_tsv(cor_annot_df, "./human_ebv_correlations.tsv")
#write_tsv(cor_annot_df, "./human_ebv_correlations_tpm10.tsv")
write_tsv(cor_annot_df, "./human_ebv_correlations_tpm5in75.tsv")

# correlation with copy number

ebv_copies <- readxl::read_excel("./ebv_copynumbers_pone.0179446.xlsx")

copies_df <- inner_join(human, ebv_copies, by = c("sampleid" = "samples")) %>%
    group_by(id, gid) %>%
    summarise(corr = broom::tidy(cor.test(normexp, `EBV load`, method = "pearson", exact = FALSE))) %>%
    ungroup() %>%
    unnest(cols = c(corr)) %>%
    select(id, gid, estimate, p = p.value)

write_tsv(copies_df, "./human_expression_ebvcopies_corr.tsv")







