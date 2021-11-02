library(SNPRelate)
library(SeqArray)

vcf_file <- commandArgs(TRUE)[1]
gds_file <- commandArgs(TRUE)[2]

# Convert VCF to GDS
seqVCF2GDS(vcf_file, gds_file, fmt.import = "GT", verbose = FALSE)

# Open GDS
gds <- seqOpen(gds_file)

# Prunning by r2 = 0.1
set.seed(1)
pruned <- snpgdsLDpruning(gds, method = "corr", ld.threshold = sqrt(0.1))

pruned_snps <- unlist(pruned, use.names = FALSE)

# Save pruned GDS
seqSetFilter(gds, variant.id = pruned_snps)

gds_out <- sub("\\.gds$", "_pruned.gds", gds_file)
seqExport(gds, gds_out)

# Delete full GDS
#unlink(gds_file)
