library(tidyverse)

list.files("./sle_variants")

sle_hg38 <- read_tsv("./sle_variants/sle_variants_hg38.bed", col_names = FALSE)

sle <- read_tsv("./sle_variants/sle_variants.tsv") %>%
    add_column(pos_hg38 = sle_hg38$X2)

chr_keys <- read_tsv("./sle_variants/dbSNP_conversion_key.txt")

dbsnp <- read_tsv("./sle_variants/sle_variants_hg38_dbSNP.vcf", comment = "##") %>%
    mutate(type = str_extract(INFO, "(?<=VC=)[^;]+")) %>%
    left_join(chr_keys, by = c("#CHROM" = "id")) %>%
    select(chr, pos = 2, snp_id = 3, ref = 4, alt = 5, type) %>%
    inner_join(sle, by = c("chr", "pos" = "pos_hg38"))

vcf <- read_tsv("./sle_variants/sle.MGB.vcf", comment = "##")
vcf_info <- select(vcf, chr = 1, pos = 2, snp_id = 3, ref = 4, alt = 5)

vcf_info %>% count(chr, pos, sort = TRUE)

sle %>% filter(chr == "chr9", pos_hg38 == 98020413)
dbsnp %>% filter(chr == "chr9", pos == 98020413)
vcf_info %>% filter(chr == "chr9", pos == 98020413)

# when computing scores, use distinct() to keep only the SNV variant
