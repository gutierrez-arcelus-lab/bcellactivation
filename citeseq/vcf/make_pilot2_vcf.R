library(IRanges)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
    
# Functions
read_geno <- function(f) {
    read_tsv(f, comment = "#", col_names = c("rsid", "chr", "pos", "gt")) |>
    filter(chr %in% c(1:22, "X")) |>
    mutate(a1 = sub("^(.)(.)$", "\\1", gt),
	   a2 = sub("^(.)(.)$", "\\2", gt)) |>
    select(-gt)
}

# read genotype data
citeseq_dir <- "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/CITEseq_pilot_2"

geno_1 <- 
    file.path(citeseq_dir, "genome_Maria_Gutierrez_v5_Full_20220718154317.txt") |>
    read_geno()

geno_2 <- 
    file.path(citeseq_dir, "genome_Roxane_Darbousset_v5_Full_20220509124117.txt") |>
    read_geno()

# merge individuals
# remove missing
# remove indels (alleles 'D'/'I')
geno_merge <- 
    inner_join(geno_1, geno_2, 
	       by = c("rsid", "chr", "pos"),
	       suffix = c("_SAMPLE1", "_SAMPLE2")) |>
    filter(grepl("^rs", rsid)) |>
    pivot_longer(-(1:3), names_to = "hap", values_to = "allele") |>
    group_by(rsid, chr, pos) |>
    filter(all(allele %in% c("A", "T", "C", "G"))) |>
    ungroup() |>
    pivot_wider(names_from = hap, values_from = allele)

# liftover positions to GRCh38
geno_bed <- geno_merge |>
    mutate(chr = paste0("chr", chr),
	   start = pos - 1L) |>
    select(chr, start, end = pos, rsid) |>
    arrange(chr, start)

bed_file_19 <- "./data/pos_23andme_hg19.bed"
bed_file_38 <- "./data/pos_23andme_hg38.bed"
chain <- "/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz"
fail <- "./data/failToLift.txt"

write_tsv(geno_bed, bed_file_19, col_names = FALSE)

command <- sprintf("liftOver %s %s %s %s", bed_file_19, chain, bed_file_38, fail)
system(command)

geno_bed38 <- read_tsv(bed_file_38, col_names = c("chr", "start", "end", "rsid"))

# Update rsid
geno_bed38 |>
    mutate(rsid = sub("^rs", "", rsid)) |>
    pull(rsid) |>
    write_lines("./data/rsid_23andMe.txt")

system("./update_rsid_23andme.sh")

updated <- read_tsv("./data/rsid_23andMe_updated.txt", col_names = c("status", "rsid_new")) |>
    mutate(rsid_new = paste0("rs", rsid_new))

geno_bed38_upd <- bind_cols(geno_bed38, updated) |>
    select(chr, start, end, rsid, rsid_new)

geno_merged_lifted <- geno_merge |>
    mutate(chr = paste0("chr", chr)) |>
    inner_join(geno_bed38_upd, by = c("chr", "rsid")) |>
    select(chr, pos = end, rsid = rsid_new, a1_SAMPLE1:a2_SAMPLE2) |>
    distinct()

# Select for 2 alleles seen
# Transform major allele to 'REF' and minor allele to 'ALT'
geno_status <- geno_merged_lifted |>
    pivot_longer(a1_SAMPLE1:a2_SAMPLE2, names_to = "h", values_to = "allele") |>
    count(chr, pos, rsid, allele) |>
    group_by(chr, pos, rsid) |>
    filter(n() == 2) |>
    arrange(chr, pos, rsid, desc(n)) |>
    mutate(status = c("REF", "ALT")) |>
    ungroup() |>
    select(-n) |>
    pivot_wider(names_from = status, values_from = allele)

# Make VCF
vcf <- geno_merged_lifted |>
    inner_join(geno_status, by = c("chr", "pos", "rsid")) |> 
    pivot_longer(a1_SAMPLE1:a2_SAMPLE2, names_to = "h", values_to = "allele") |>
    mutate(ds = case_when(allele == REF ~ 0, 
			  allele == ALT ~ 1)) |>
    select(-allele) |>
    group_by(chr, pos, rsid) |>
    mutate(ac = sum(ds),
	   an = n()) |>
    ungroup() |>
    pivot_wider(names_from = h, values_from = ds) |>
    mutate(af = ac/an,
	   info = paste0("AC=", ac, ";AN=", an, ";AF=", af, ";VT=SNP;NS=2"),
	   QUAL = ".", 
	   FILTER = "PASS", 
	   FORMAT = "GT",
	   SAMPLE1 = paste(pmin(a1_SAMPLE1, a2_SAMPLE1), pmax(a1_SAMPLE1, a2_SAMPLE1), sep = "/"),
	   SAMPLE2 = paste(pmin(a1_SAMPLE2, a2_SAMPLE2), pmax(a1_SAMPLE2, a2_SAMPLE2), sep = "/")) |>
    select(`#CHROM` = chr, POS = pos, ID = rsid, REF, ALT,
	   QUAL, FILTER, INFO = info, FORMAT, SAMPLE1, SAMPLE2) |>
    arrange(`#CHROM`, POS)


out_all <- "./data/pilot2_genotypes.vcf"
file.create(out_all)
write_lines("##fileformat=VCFv4.1", out_all)
write_lines('##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">',
	    out_all, append = TRUE)
write_lines('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
	    out_all, append = TRUE)
write_lines('##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">',
	    out_all, append = TRUE)
write_lines('##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant">',
	    out_all, append = TRUE)
write_lines('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">',
	    out_all, append = TRUE)
write_lines('##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased Genotype">',
	    out_all, append = TRUE)
write_tsv(vcf, out_all, append = TRUE, col_names = TRUE)
unlink(paste0(out_all, ".gz"))
system(paste("bgzip", out_all))
system(paste0("tabix -p vcf ", out_all, ".gz"))


# select SNPs within gene coordinates
# (not implemented in final version of VCF)
annot <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf") |>
    read_tsv(comment = "#", col_names = FALSE)

annot_genes <- annot |>
    filter(X1 %in% paste0("chr", c(1:22, "X")),
	   X3 == "gene",
	   grepl("gene_type\\s\"(protein_coding|lncRNA)\"", X9)) |>
    distinct(chr = X1, start = X4, end = X5)

vcf_genes <- vcf |>
    inner_join(annot_genes, 
	       join_by(`#CHROM` == chr, between(POS, start, end)))
