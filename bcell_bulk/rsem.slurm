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

# Genome sequence
GENOME=$HOME/lupus/data/GRCh38.primary_assembly.genome.fa

# Gencode annotations
GTF=$HOME/lupus/data/gencode.v38.primary_assembly.annotation.gtf

# Results
OUT=$HOME/lupus/data/gencode

# Temporarily unzip because rsem does not take zipped files
zcat $GENOME.gz > $GENOME
zcat $GTF.gz > $GTF

# Run rsem to make our transcript fasta;
# given the genome sequence and annotations
rsem-prepare-reference --gtf $GTF $GENOME $OUT

# Remove unnecessary output files
rm $GENOME $GTF
rm $OUT.grp $OUT.ti $OUT.chrlist $OUT.seq $OUT.idx.fa $OUT.n2g.idx.fa
