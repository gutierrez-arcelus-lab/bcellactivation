library(tidyverse)

report <- 
    file.path("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405",
	      "GCF_000001405.39_GRCh38.p13",
	      "GCF_000001405.39_GRCh38.p13_assembly_report.txt") %>%
    read_tsv(comment = "#", col_names = FALSE) %>%
    filter(X4 == "Chromosome") %>%
    filter(X1 %in% c(1:22, "X")) %>%
    select(chr = X10, id = X7)

sle <-
    "./sle_variants/sle_variants_hg38.txt" %>%
    read_tsv(col_names = c("chr", "pos")) %>%
    left_join(report, by = "chr") %>%
    select(chr = id, pos)

write_tsv(report, "./sle_variants/dbSNP_conversion_key.txt")

write_tsv(sle, 
	  "./sle_variants/sle_variants_hg38_dbSNP.txt",
	  col_names = FALSE)
