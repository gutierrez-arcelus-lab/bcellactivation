#!/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=2:00:00
#SBATCH --mem=48gb
#SBATCH -p bch-compute
#SBATCH --job-name=STAR-index
#SBATCH --mail-user=vitor.RezendeDaCostaAguiar@childrens.harvard.edu
#SBACTH --mail-type=END,FAIL
#SBATCH -o STARindex.o%j

source /programs/biogrids.shrc
export STAR_X=2.7.9a

GENOME=/reference_databases/ReferenceGenome/STAR/GRCh38.v32/GRCh38.primary_assembly.genome.fa
GTF=/reference_databases/ReferenceGenome/STAR/GRCh38.v32/gencode.v32.primary_assembly.annotation.gtf 
OUT=/lab-share/IM-Gutierrez-e2/Public/vitor/starindex

mkdir -p $OUT

STAR --runThreadN $SLURM_CPUS_PER_TASK \
    --runMode genomeGenerate \
    --genomeDir $OUT \
    --genomeFastaFiles $GENOME \
    --sjdbGTFfile $GTF
