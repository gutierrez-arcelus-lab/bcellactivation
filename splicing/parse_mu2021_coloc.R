library(tidyverse)
library(readxl)

# intron-gene map
intron_gene <- read_tsv("./mu2021_data/sqtl_all.tsv", col_types = c(sid = "c")) %>%
    select(intron = pid, gene = genes) %>%
    mutate(intron = str_replace(intron, "_(?=[+-])", ":"),
	   intron = str_remove(intron, "^chr"))

# Co-localization results
coloc_file <- "./mu2021_data/colocalization_results.xlsx"

datasets <- excel_sheets(coloc_file)[-7]

coloc_eqtl <- grep("eQTL", datasets, value = TRUE) %>%
    map_df(~read_excel(coloc_file, ., skip = 1)) %>%
    filter(trait == "SLE") %>%
    separate(colocSnp, c("chr_gwas", "pos_gwas"), sep = ":", convert = TRUE) %>%
    select(study, cell, gene, chr_gwas, pos_gwas, coloc_locus = colocLoci, pp4 = PP.H4.abf) %>%
    group_by(study, cell, gene) %>%
    slice(which.max(pp4)) %>%
    ungroup()

coloc_sqtl <- grep("sQTL", datasets, value = TRUE) %>%
    map_df(~read_excel(coloc_file, ., skip = 1)) %>%
    filter(trait == "SLE") %>%
    separate(colocSnp, c("chr_gwas", "pos_gwas"), sep = ":", convert = TRUE) %>%
    select(study, cell, intron, chr_gwas, pos_gwas, coloc_locus = colocLoci, pp4 = PP.H4.abf) %>%
    left_join(intron_gene, by = "intron") %>%
    filter(!is.na(gene)) %>%
    group_by(study, cell, gene) %>%
    slice(which.max(pp4)) %>%
    ungroup()


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
	   p = `P-VALUE`) %>%
    mutate(snp_current = ifelse(!is.na(snp_current), paste0("rs", snp_current), snp_current),
	   pos = as.integer(pos)) %>%
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
    select(chr, pos = pos_hg19, reported_gene, snp, p) %>%
    mutate(chr = ifelse(is.na(chr) & reported_gene == "BANK1", sub("^([^:]+):.*$", "\\1", snp), chr),
	   pos = ifelse(is.na(pos) & reported_gene == "BANK1", sub("^[^:]+:(.+)$", "\\1", snp), pos),
	   pos = as.integer(pos))

write_tsv(gwas_hg19, "./plot_data/sle_gwas.tsv")
write_tsv(coloc_eqtl, "./plot_data/coloc_eqtl.tsv")
write_tsv(coloc_sqtl, "./plot_data/coloc_sqtl.tsv")

