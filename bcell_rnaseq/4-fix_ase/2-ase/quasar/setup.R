library(tidyverse)

if (!file.exists("data")) dir.create("data")

# Make a VCF with variants originally selected for ASE
# to compute their AFs in MGB
vcf <-
    "../../0-genotypes/data/allchr.mgb.vcf.gz" |>
    read_tsv(comment = "##")

variant_vcf <- 
    vcf |>
    select(`#CHROM`:INFO) |>
    mutate(FILTER = ".", 
	   INFO = ".") 
   
write_lines("##fileformat=VCFv4.1", "./data/variants.vcf")
write_tsv(variant_vcf, "./data/variants.vcf", append = TRUE, col_names = TRUE)
system("bgzip ./data/variants.vcf")
system("tabix -p vcf ./data/variants.vcf.gz")