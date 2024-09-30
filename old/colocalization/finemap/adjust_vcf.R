library(tidyverse)

cmdargs <- commandArgs(TRUE)
vcfin <- cmdargs[1]
var_file <- cmdargs[2]
vcfout <- cmdargs[3]

Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)

## At multiallelic variants, remove allele with 0 counts in this sample
## and recode alleles so some will be now biallelic.
## Then we don't lose these variants when we 
variants <- read_lines(var_file)

tmp <- read_lines(vcfin, n_max = 5000)
header <- keep(tmp, grepl("^##", tmp))

# deactivate adjustment because need to deal of alleles encoded as 2.
# Need to transform them to "1" after filters
adj_vcf <- read_tsv(vcfin, comment = "##") |>
    #extract(INFO, "ac", "^AC=([^;]+)", remove = FALSE) |>
    #separate_rows(c(ALT, ac), sep = ",", convert = TRUE) |>
    #filter(ac != 0) |>
    #select(-ac) |>
    #group_by(across(c(-ALT))) |>
    #summarise(ALT = paste(ALT, collapse = ",")) |>
    #ungroup() |>
    mutate(QUAL = ".", FILTER = ".", INFO = ".",
	   ID = paste(`#CHROM`, POS, REF, ALT, sep = ":")) |>
    select(`#CHROM`, POS, ID, REF, ALT, everything()) |>
    filter(ID %in% variants)

vcfout_tmp <- sub("\\.gz$", "", vcfout)
write_lines(header, vcfout_tmp)
write_tsv(adj_vcf, vcfout_tmp, col_names = TRUE, append = TRUE)
unlink(vcfout)
#unlink(paste0(vcfout, ".tbi"))
system(paste("bgzip", vcfout_tmp))
#system(paste("tabix -p vcf", vcfout))
