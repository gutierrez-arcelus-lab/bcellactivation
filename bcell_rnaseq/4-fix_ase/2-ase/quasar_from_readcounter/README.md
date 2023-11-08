After running QuASAR, I noticed that it throws away any variants with a 'N' in the CIGAR string.
That means discarding, for example, exonic variants that are spliced out in some alignments.

Therefore, we decided to run QuASAR using data from GATK ASEReadCounter, thus not using
the de novo genotyping from QuASAR.
