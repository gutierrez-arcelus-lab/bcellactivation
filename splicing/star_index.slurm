#!/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=2:00:00
#SBATCH --mem=36gb
#SBATCH -p bch-compute-pe
#SBATCH --job-name=STARidx
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBACTH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/STAR-index.o%j

source /programs/biogrids.shrc
export STAR_X=2.7.9a

LABSHR=/lab-share/IM-Gutierrez-e2/Public 
GENOME=${LABSHR}/References/Genomes/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa
GTF=${LABSHR}/References/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf
OUT=${SLURM_SUBMIT_DIR}/starindex_51bp

mkdir -p $OUT

STAR --runMode genomeGenerate \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $OUT \
    --genomeFastaFiles $GENOME \
    --sjdbGTFfile $GTF \
    --sjdbOverhang 50
