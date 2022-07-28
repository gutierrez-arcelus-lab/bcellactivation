library(IRanges)
library(tidyverse)
    
# Functions
read_geno <- function(f) {
    read_tsv(f, comment = "#", col_names = c("rsid", "chr", "pos", "gt")) %>%
    filter(chr %in% c(1:22, "X")) %>%
    mutate(a1 = sub("^(.)(.)$", "\\1", gt),
	   a2 = sub("^(.)(.)$", "\\2", gt)) %>%
    select(-gt)
}

# read genotype data
citeseq_dir <- "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/CITEseq_pilot_2"

geno_1 <- 
    file.path(citeseq_dir, "genome_Maria_Gutierrez_v5_Full_20220718154317.txt") %>%
    read_geno()

geno_2 <- 
    file.path(citeseq_dir, "genome_Roxane_Darbousset_v5_Full_20220509124117.txt") %>%
    read_geno()

# merge individuals
# remove missing
# select SNPs with up to 2 different alleles
geno_merge <- 
    inner_join(geno_1, geno_2, 
	       by = c("rsid", "chr", "pos"),
	       suffix = c("_SAMPLE1", "_SAMPLE2")) %>%
    pivot_longer(-(1:3), names_to = "hap", values_to = "allele") %>%
    group_by(rsid, chr, pos) %>%
    filter(!any(allele == "-")) %>%
    filter(n_distinct(allele) <= 2) %>%
    ungroup() %>%
    pivot_wider(names_from = hap, values_from = allele)

db_chrs_ids_19 <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.25_GRCh37.p13_assembly_report.txt" %>%
    read_tsv(comment = "#", col_names = FALSE) %>%
    filter(X1 %in% c(1:22, "X")) %>%
    distinct(X1, X7)

db_chrs_ids_38 <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.39_GRCh38.p13_assembly_report.txt" %>%
    read_tsv(comment = "#", col_names = FALSE) %>%
    filter(X1 %in% c(1:22, "X")) %>%
    distinct(X1, X7)

geno_merge %>%
    left_join(db_chrs_ids_19, by = c("chr" = "X1")) %>%
    select(chr = X7, pos) %>%
    arrange(chr) %>%
    write_tsv("./pos_23andme.txt", col_names = FALSE)

# liftover positions to GRCh38
geno_bed <- geno_merge %>%
    mutate(chr = paste0("chr", chr),
	   start = pos - 1L) %>%
    select(chr, start, end = pos, rsid) %>%
    arrange(chr, start)

bed_file_19 <- "./pos_23andme_hg19.bed"
bed_file_38 <- "./pos_23andme_hg38.bed"
chain <- "/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz"
fail <- "./failToLift.txt"

write_tsv(geno_bed, bed_file_19, col_names = FALSE)

command <- sprintf("liftOver %s %s %s %s", bed_file_19, chain, bed_file_38, fail)
system(command)

geno_bed38 <- read_tsv(bed_file_38, col_names = c("chr", "start", "end", "rsid"))

geno_bed38 %>%
    mutate(chr = sub("chr", "", chr)) %>%
    left_join(db_chrs_ids_38, by = c("chr" = "X1")) %>%
    select(chr = X7, pos = end) %>%
    arrange(chr) %>%
    write_tsv("./pos_23andme_hg38.txt", col_names = FALSE)

# Extract positions from dbSNP on assemblies hg19 and GRCh38
# Run:
# ./subset_dbsnp.sh

dbsnp_19 <- "./dbsnp_23andme.vcf" %>%
    read_tsv(comment = "##", col_select = c(1:5, 8)) %>%
    filter(grepl("VC=SNV", INFO)) %>%
    left_join(db_chrs_ids_19, by = c("#CHROM" = "X7")) %>%
    select(`#CHROM` = X1, POS, ID, REF, ALT)

dbsnp_38 <- "./dbsnp_23andme_hg38.vcf" %>%
    read_tsv(comment = "##", col_select = c(1:5, 8)) %>%
    filter(grepl("VC=SNV", INFO)) %>%
    left_join(db_chrs_ids_38, by = c("#CHROM" = "X7")) %>%
    select(`#CHROM` = X1, POS, ID, REF, ALT)

dbsnp <- 
    inner_join(dbsnp_19, dbsnp_38, 
	       by = c("#CHROM", "ID"),
	       suffix = c("_19", "_38")) %>%
    filter(REF_19 == REF_38, ALT_19 == ALT_38) %>%
    select(`#CHROM`, POS_19, POS_38, ID, REF = REF_19, ALT = ALT_19)


# Select the ALT allele present in the data when multiple ALT are possible
variants_df <- geno_merge %>%
    pivot_longer(-(1:3), names_to = "hap", values_to = "allele") %>%
    distinct(rsid, chr, pos, allele) %>%
    group_by(rsid, chr, pos) %>%
    summarise(alleles = list(allele)) %>%
    ungroup() %>%
    inner_join(dbsnp, by = c("chr" = "#CHROM", "pos" = "POS_19", "rsid" = "ID")) %>%
    separate_rows(ALT, sep = ",") %>%
    mutate(i = map2_chr(ALT, alleles, ~ifelse(.x %in% unlist(.y), 1L, 0L))) %>%
    group_by(rsid, chr, pos) %>%
    filter((all(i == 0L) & row_number() == 1L) | (any(i == 1L) & i == 1L)) %>%
    ungroup() %>%
    select(chr, pos, POS_38, rsid, REF, ALT) %>%
    arrange(chr, pos)

# Make VCF
vcf <- geno_merge %>%
    inner_join(variants_df, by = c("chr", "pos", "rsid")) %>%
    mutate_at(vars(a1_SAMPLE1:a2_SAMPLE2), 
	      ~case_when(. == REF ~ 0L, . == ALT ~ 1L, TRUE ~ NA_integer_)) %>%
    pivot_longer(a1_SAMPLE1:a2_SAMPLE2, names_to = "hap", values_to = "allele") %>%
    group_by(rsid, chr, pos, REF, ALT) %>%
    mutate(ac = sum(allele)) %>%
    ungroup() %>%
    pivot_wider(names_from = hap, values_from = allele) %>%
    mutate(an = 4L, 
	   af = ac/an,
	   info = paste0("AC=", ac, ";AN=", an, ";AF=", af, ";VT=SNP;NS=2"),
	   QUAL = ".", 
	   FILTER = "PASS", 
	   FORMAT = "GT",
	   SAMPLE1 = paste(pmin(a1_SAMPLE1, a2_SAMPLE1), pmax(a1_SAMPLE1, a2_SAMPLE1), sep = "/"),
	   SAMPLE2 = paste(pmin(a1_SAMPLE2, a2_SAMPLE2), pmax(a1_SAMPLE2, a2_SAMPLE2), sep = "/")) %>%
    select(`#CHROM` = chr, POS = POS_38, ID = rsid, REF, ALT,
	   QUAL, FILTER, INFO = info, FORMAT, SAMPLE1, SAMPLE2) %>%
    mutate(`#CHROM` = paste0("chr", `#CHROM`)) %>%
    arrange(`#CHROM`, POS)


out_all <- "./demuxlet/samples_allsnps.vcf"
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


# select exonic SNPs
annot <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf") %>%
    read_tsv(comment = "#", col_names = FALSE)

annot_exons <- annot %>%
    filter(X1 %in% paste0("chr", c(1:22, "X")),
	   X3 == "exon",
	   grepl("gene_type\\s\"(protein_coding|lncRNA)\"", X9)) %>%
    distinct(chr = X1, start = X4, end = X5) %>%
    mutate(chr = sub("chr", "", chr))

exons_ranges <- annot_exons %>%
    split(.$chr) %>%
    map(~IRanges(.$start, .$end)) %>%
    .[c(1:22, "X")]

pos_list <- vcf %>%
    split(.$`#CHROM`) %>%
    map(~IRanges(.$POS, .$POS)) %>%
    .[c(1:22, "X")]

snps_overlap <- map2(pos_list, exons_ranges, subsetByOverlaps) %>%
    map(as.data.frame) %>%
    bind_rows(.id = "chr") %>%
    as_tibble() %>%
    select(-width)

vcf_exons <- snps_overlap %>%
    select(chr, pos = start) %>%
    inner_join(vcf38_filt, ., by = c("#CHROM" = "chr", "POS" = "pos")) %>%
    mutate(`#CHROM` = paste0("chr", `#CHROM`)) %>%
    arrange(`#CHROM`)

out <- "./samples.vcf"
file.create(out)
write_lines("##fileformat=VCFv4.1", out)
write_tsv(vcf_exons, out, append = TRUE, col_names = TRUE)
unlink("./samples.vcf.gz")
system(paste("bgzip", out))
system(paste0("tabix -p vcf ", out, ".gz"))


