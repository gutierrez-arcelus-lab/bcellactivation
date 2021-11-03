library(SNPRelate)
library(SeqArray)

vcf_file <- commandArgs(TRUE)[1]
gds_file <- commandArgs(TRUE)[2]

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

gds_out <- sub("\\.gds$", ".pruned.gds", gds_file)
seqExport(gds, gds_out)

vcf_out <- sub("gds$", "vcf", gds_out)
seqGDS2VCF(gds_out, vcf_out)
