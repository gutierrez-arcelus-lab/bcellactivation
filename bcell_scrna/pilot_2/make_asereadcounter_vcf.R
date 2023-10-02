library(tidyverse)

vcf_file <- "./demuxlet/samples_allsnps.vcf.gz"

vcf_header <- read_lines(vcf_file) %>%
    keep(~grepl("^##", .))

vcf <- read_tsv(vcf_file, comment = "##") %>%
    filter(SAMPLE1 != SAMPLE2) %>%
    select(-SAMPLE2) %>%
    mutate(SAMPLE1 = "0/1")

out_file <- "./demuxlet/samples_allsnps_asereadc.vcf"
file.create(out_file)
write_lines(vcf_header, out_file, append = TRUE)
write_tsv(vcf, out_file, col_names = TRUE, append = TRUE)
unlink(paste0(out_file, ".gz"))
system(paste("bgzip", out_file))
system(paste0("tabix -f -p vcf ", out_file, ".gz"))
