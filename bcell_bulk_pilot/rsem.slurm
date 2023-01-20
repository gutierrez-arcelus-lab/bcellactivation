#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=RSEM
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/RSEM.%j

source /programs/biogrids.shrc
export RSEM_X=1.2.25

LABSHR=/lab-share/IM-Gutierrez-e2/Public/References

# Genome sequence
GENOME=${LABSHR}/Genomes/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa

# Gencode annotations
GTF=${LABSHR}/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf

# Results
OUT=${LABSHR}/Annotations/hsapiens/gencode.v38.primary_assembly

# Run rsem to make our transcript fasta;
# given the genome sequence and annotations
rsem-prepare-reference --gtf $GTF $GENOME $OUT

# Remove unnecessary output files
rm $OUT.grp $OUT.ti $OUT.chrlist $OUT.seq $OUT.idx.fa $OUT.n2g.idx.fa
