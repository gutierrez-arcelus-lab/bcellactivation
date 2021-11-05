library(SNPRelate)
library(SeqArray)

vcf_file <- commandArgs(TRUE)[1]
out_prefix <- commandArgs(TRUE)[2]
gds_file <- paste0(out_prefix, ".gds")

# Convert VCF to GDS
seqVCF2GDS(vcf_file, gds_file, fmt.import = "GT", verbose = FALSE)

# Open GDS
gds <- seqOpen(gds_file)

# Prunning by r2 = 0.2
set.seed(1)
pruned <- snpgdsLDpruning(gds,
			  slide.max.bp = 5e4,
			  method = "corr", 
			  ld.threshold = sqrt(0.2))

pruned_snps <- unlist(pruned, use.names = FALSE)

# Save pruned GDS
seqSetFilter(gds, variant.id = pruned_snps)

gds_out <- paste0(out_prefix, ".pruned.gds")
seqExport(gds, gds_out)

vcf_out <- paste0(out_prefix, ".vcf.gz")
seqGDS2VCF(gds_out, vcf_out)
