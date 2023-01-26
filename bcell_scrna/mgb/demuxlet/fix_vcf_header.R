vcf <- "./data/allchr.mgb.sorted.vcf.gz" 

header <- system(sprintf("bcftools view -h %s", vcf), intern = TRUE)

contigs <- grep("^##contig", header)

sorted_contigs <- sort(header[contigs])

tmp_header <- header[-contigs]

final_header <- c(tmp_header[-length(tmp_header)], sorted_contigs, tmp_header[length(tmp_header)])

writeLines(final_header, "./data/header.txt")

