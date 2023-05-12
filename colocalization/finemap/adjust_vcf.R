library(tidyverse)

cmdargs <- commandArgs(TRUE)
vcfin <- cmdargs[1]
vcfout <- cmdargs[2]

Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)

## At multiallelic variants, remove allele with 0 counts in this sample
## and recode alleles so some will be now biallelic
tmp <- read_lines(vcfin, n_max = 5000)
header <- keep(tmp, grepl("^##", tmp))

adj_vcf <- read_tsv(vcfin, comment = "##") |>
    extract(INFO, "ac", "^AC=([^;]+)", remove = FALSE) |>
    separate_rows(c(ALT, ac), sep = ",") |>
    filter(ac != "0") |>
    select(-ac) |>
    group_by(across(c(-ALT))) |>
    summarise(ALT = paste(ALT, collapse = ",")) |>
    ungroup() |>
    mutate(QUAL = ".", FILTER = ".", INFO = ".",
	   ID = paste(`#CHROM`, POS, REF, ALT, sep = ":")) |>
    select(`#CHROM`, POS, ID, REF, ALT, everything())

write_lines(header, vcfout)
write_tsv(adj_vcf, vcfout, col_names = TRUE, append = TRUE)
unlink(paste0(vcfout, ".gz"))
system(paste("bgzip", vcfout))
system(paste0("tabix -p vcf ", vcfout, ".gz"))

