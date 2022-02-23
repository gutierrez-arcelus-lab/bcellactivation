library(tidyverse)
library(Biostrings)

human <- readDNAStringSet("./gencode.transcripts.fa")
ebv <- readDNAStringSet("./ebv_genbank.fasta")

locus <- str_extract(names(ebv), "(?<=locus_tag=)[^]]+")
seqid <- sub("^lcl\\|([^ ]+).+$", "\\1", names(ebv))

names(ebv) <- paste(locus, seqid, sep = ":")

human_ebv <- c(human, ebv)

writeXStringSet(human_ebv, "./gencode_plusEBV.transcripts.fa")
