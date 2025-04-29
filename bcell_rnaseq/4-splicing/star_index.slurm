#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=48gb
#SBATCH -p bch-compute
#SBATCH --job-name=STARidx
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBACTH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/STAR-index.o%A-%a.out
#SBATCH -e /temp_work/ch229163/log/STAR-index.o%A-%a.err

source /programs/biogrids.shrc
export STAR_X=2.7.11a

LABSHR=/lab-share/IM-Gutierrez-e2/Public 
GENOME=${LABSHR}/References/Genomes/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa
GTF=${LABSHR}/References/Annotations/hsapiens/gencode.v41.primary_assembly.annotation.gtf
OH=149

OUT=${SLURM_SUBMIT_DIR}/index/STAR

mkdir -p $OUT

STAR --runMode genomeGenerate \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $OUT \
    --genomeFastaFiles $GENOME \
    --sjdbGTFfile $GTF \
    --sjdbOverhang $OH
