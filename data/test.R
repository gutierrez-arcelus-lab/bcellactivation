library(tidyverse)

dbsnp <- read_tsv("./dbsnp_155.hg38pri.vcf.gz", comment = "##", col_select = 1:3)
kgp <- read_tsv("./1000g.gatkbundle.vcf.gz", comment = "##", col_select = 1:3)
indel <- read_tsv("./knownindels_gatkbundle.hg38pri.vcf.gz", comment = "##", col_select = 1:3)
mills <- read_tsv("./mills_1000g_indels_gatkbundle.hg38pri.vcf.gz", comment = "##", col_select = 1:3)

library(Biostrings)

fasta <- readDNAStringSet("./GRCh38.primary_assembly.genome.fa")

names(fasta) <- sub(" (.+)$", "", names(fasta))

writeXStringSet(fasta, "./GRCh38.primary_assembly.genome.fixnames.fa")
