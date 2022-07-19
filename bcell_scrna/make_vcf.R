library(tidyverse)
    
# Functions
read_geno <- function(f)
    read_tsv(f, comment = "#", col_names = c("rsid", "chr", "pos", "gt")) %>%
    filter(chr %in% c(1:22, "X")) %>%
    mutate(a1 = sub("^(.)(.)$", "\\1", gt),
	   a2 = sub("^(.)(.)$", "\\2", gt)) %>%
    select(-gt)


citeseq_dir <- "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/CITEseq_pilot_2"

geno_1 <- 
    file.path(citeseq_dir, "genome_Maria_Gutierrez_v5_Full_20220718154317.txt") %>%
    read_geno()

geno_2 <- 
    file.path(citeseq_dir, "genome_Roxane_Darbousset_v5_Full_20220509124117.txt") %>%
    read_geno()

geno_merge <- inner_join(geno_1, geno_2, by = c("rsid", "chr", "pos")) %>%
    filter(a1.x != "-", a2.x != "-", a1.y != "-", a2.y != "-")
   
geno_merge %>%
    pull(rsid) %>%
    write_lines(file.path(citeseq_dir, "./snps_23andMe.txt"))

# run:
# ./subset_dbsnp.sh

db_report <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.25_GRCh37.p13_assembly_report.txt" %>%
    read_tsv(comment = "#", col_names = FALSE) %>%
    filter(X1 %in% c(1:22, "X")) %>%
    distinct(X1, X7)

dbsnp <- file.path(citeseq_dir, "dbsnp.vcf") %>%
    read_tsv(comment = "##", col_select = 1:5) %>%
    filter(`#CHROM` %in% db_report$X7) %>%
    select(ID, REF, ALT) %>%
    filter(REF %in% c("A", "G", "C", "T"), 
	   ALT %in% c("A", "G", "C", "T"))

vcf <- inner_join(geno_merge, dbsnp, by = c("rsid" = "ID")) %>%
    mutate_at(vars(a1.x:a2.y), ~case_when(. == REF ~ 0L, . == ALT ~ 1L, TRUE ~ NA_integer_)) %>%
    drop_na() %>%
    unite("sample_1", c("a1.x", "a2.x"), sep = "/") %>%
    unite("sample_2", c("a1.y", "a2.y"), sep = "/") %>%
    filter(sample_1 != sample_2) %>%
    mutate(QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT") %>%
    select(`#CHROM` = chr, POS = pos, ID = rsid, REF, ALT,
	   QUAL, FILTER, INFO, FORMAT,
	   sample_1, sample_2)


# select exonic SNPs
library(IRanges)

annot <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References", 
	      "Annotations/hsapiens/gencode.v19.chr_patch_hapl_scaff.annotation.gtf") %>%
    read_tsv(comment = "#", col_names = FALSE)

annot_exons <- annot %>%
    filter(X1 %in% paste0("chr", c(1:22, "X")),
	   X3 == "exon",
	   grepl("gene_type\\s\"(protein_coding|lincRNA)\"", X9)) %>%
    distinct(chr = X1, start = X4, end = X5) %>%
    mutate(chr = sub("chr", "", chr))

exons_ranges <- annot_exons %>%
    split(.$chr) %>%
    map(~IRanges(.$start, .$end)) %>%
    .[c(1:22, "X")]

pos_list <- vcf %>%
    split(.$`#CHROM`) %>%
    map(~IRanges(.$POS, .$POS)) %>%
    .[c(1:22, "X")]

vcf_exons <- map2(pos_list, exons_ranges, subsetByOverlaps) %>%
    map(as.data.frame) %>%
    bind_rows(.id = "chr") %>%
    as_tibble() %>%
    select(chr, pos = start) %>%
    inner_join(vcf, ., by = c("#CHROM" = "chr", "POS" = "pos"))


# Subsample 
vcf_exons %>%
    sample_n(100) %>%
    mutate(`#CHROM` = factor(`#CHROM`, levels = str_sort(unique(`#CHROM`), numeric = TRUE))) %>%
    arrange(`#CHROM`, POS)
