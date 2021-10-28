#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=RSEM
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /home/ch229163/ase/log/RSEM.%j

source /programs/biogrids.shrc
export RSEM_X=1.2.25

GENOME=$HOME/ase/data/GRCh38.primary_assembly.genome.fa
GTF=$HOME/ase/data/gencode.v38.primary_assembly.annotation.gtf
OUT=$HOME/ase/data/gencode

zcat $GENOME.gz > $GENOME
zcat $GTF.gz > $GTF

rsem-prepare-reference --gtf $GTF $GENOME $OUT

rm $GENOME $GTF
rm $OUT.grp $OUT.ti $OUT.chrlist $OUT.seq $OUT.idx.fa $OUT.n2g.idx.fa
