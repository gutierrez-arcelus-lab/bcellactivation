library(tidyverse)
library(leafcutter)
library(readxl)

exon_table <- read_tsv("./hg19_exon.txt.gz")

sqtl_dir <- "/lab-share/IM-Gutierrez-e2/Public/External_datasets/Mu2021/eqtl_sqtl_summ_stats"

qtl_df <- list.files(sqtl_dir, recursive = TRUE, full.names = TRUE) %>%
    grep("sQTL", ., value = TRUE) %>%
    setNames(., str_remove(., sqtl_dir)) %>%
    setNames(., str_remove(names(.), "^/")) %>%
    setNames(., str_remove(names(.), "\\.txt\\.gz$")) %>%
    map_df(~read_tsv(., col_types = "cidddcdddddd") %>%
	   mutate(pid = str_replace(pid, ":(?=[-+])", "_"),
		  pid = paste0("chr", pid)), 
	   .id = "source") %>%
    separate(source, c("source", "dataset"), sep = "/") %>%
    mutate(source = str_remove(source, "_sQTL"),
	   dataset = ifelse(is.na(dataset) & source == "DGN", "Whole Blood", dataset))

intron_meta <- get_intron_meta(qtl_df$pid)

clu_gene_map <- map_clusters_to_genes(intron_meta, exon_table) %>%
    as_tibble()

qtl_annot <- qtl_df %>%
    mutate(clu = sub("^(chr[^:]+):[^:]+:[^:]+:(.+)$", "\\1:\\2", pid)) %>%
    left_join(clu_gene_map, by = "clu") %>%
    select(-clu)

# intron-gene map
intron_gene <- qtl_annot %>%
    select(intron = pid, gene = genes) %>%
    mutate(intron = str_replace(intron, "_(?=[+-])", ":"),
	   intron = str_remove(intron, "^chr")) %>%
    distinct()

# Co-localization results
coloc_file <- "./colocalization_results.xlsx"

datasets <- excel_sheets(coloc_file)[-7]

coloc_eqtl <- grep("eQTL", datasets, value = TRUE) %>%
    map_df(~read_excel(coloc_file, ., skip = 1)) %>%
    filter(trait == "SLE") %>%
    separate(colocSnp, c("chr", "colocPos"), sep = ":", convert = TRUE, remove = FALSE) %>%
    select(study, cell, chr, colocSnp, colocLoci, colocPos, gene, pp4 = PP.H4.abf) %>%
    arrange(chr, colocPos) %>%
    mutate(chr = paste0("chr", chr))

coloc_sqtl <- grep("sQTL", datasets, value = TRUE) %>%
    map_df(~read_excel(coloc_file, ., skip = 1)) %>%
    filter(trait == "SLE") %>%
    separate(colocSnp, c("chr", "colocPos"), sep = ":", convert = TRUE, remove = FALSE) %>%
    select(study, cell, chr, colocSnp, colocLoci, colocPos, intron, pp4 = PP.H4.abf) %>%
    left_join(intron_gene, by = "intron") %>%
    arrange(chr, colocPos) %>%
    mutate(chr = paste0("chr", chr))

write_tsv(coloc_eqtl, "./filtered_colocs_eQTL_MuEtAl.tsv")
write_tsv(coloc_sqtl, "./filtered_colocs_sQTL_MuEtAl.tsv")



# SLE GWAS data
gwas <- "../data/gwas_catalog_v1.0.2.tsv" %>%
    read_tsv(col_types = c(.default = "c", "P-VALUE" = "d")) %>%
    filter(`FIRST AUTHOR` == "Morris DL") %>%
    select(chr = CHR_ID,
	   pos = CHR_POS, 
	   reported_gene = `REPORTED GENE(S)`,
	   mapped_gene = MAPPED_GENE,
	   snp = SNPS,
	   snp_current = SNP_ID_CURRENT,
	   context = CONTEXT,
	   population = `P-VALUE (TEXT)`,
	   p = `P-VALUE`) %>%
    mutate(snp_current = ifelse(!is.na(snp_current), paste0("rs", snp_current), snp_current),
	   pos = as.integer(pos),
	   population = str_remove_all(population, "[()]"),
	   population = ifelse(is.na(population), "Meta", population)) %>%
    arrange(as.numeric(chr), pos)

# liftOver GWAS
bed <- gwas %>%
    select(chr, pos) %>%
    mutate(chr = as.integer(chr),
	   pos = as.integer(pos),
	   end = pos + 1L) %>%
    select(chr, start = pos, end) %>%
    arrange(chr, start) %>%
    mutate(chr = paste0("chr", chr)) %>%
    filter(!is.na(start))

bed38 <- "./gwas_cat_hg38.bed"
write_tsv(bed, "./gwas_cat_hg38.bed", col_names = FALSE)

bed19 <- "./gwas_cat_hg19.bed"
chain <- "/reference_databases/ReferenceGenome/liftover_chain/hg38/hg38ToHg19.over.chain.gz"
fail <- "gwas_cat_failTolift.txt"
command <- sprintf("liftOver %s %s %s %s", bed38, chain, bed19, fail)
system(command)

map_positions <- read_tsv(bed19, col_names = c("chr", "start", "end"), col_types = c("cii")) %>%
    select(chr, pos_hg19 = start) %>%
    bind_cols(select(bed, pos_hg38 = start)) %>%
    mutate(chr = as.character(parse_number(chr)))

gwas_hg19 <- gwas %>%
    left_join(map_positions, by = c("chr", "pos" = "pos_hg38")) %>%
    distinct() %>%
    select(chr, pos = pos_hg19, reported_gene, snp, population, p) %>%
    mutate(chr = ifelse(is.na(chr) & reported_gene == "BANK1", sub("^([^:]+):.*$", "\\1", snp), chr),
	   pos = ifelse(is.na(pos) & reported_gene == "BANK1", sub("^[^:]+:(.+)$", "\\1", snp), pos),
	   pos = as.integer(pos))

write_tsv(gwas_hg19, "./plot_data/sle_gwas.tsv")
write_tsv(coloc_eqtl, "./plot_data/coloc_eqtl.tsv")
write_tsv(coloc_sqtl, "./plot_data/coloc_sqtl.tsv")

