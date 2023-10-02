#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=36gb
#SBATCH -p bch-compute
#SBATCH --job-name=STARidx
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBACTH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/STAR-index.o%A-%a.out
#SBATCH -e /temp_work/ch229163/log/STAR-index.o%A-%a.err

source /programs/biogrids.shrc
STAR=/lab-share/IM-Gutierrez-e2/Public/vitor/software/STAR/STAR_2.7.11a/Linux_x86_64_static/STAR

LABSHR=/lab-share/IM-Gutierrez-e2/Public 
GENOME=${LABSHR}/References/Genomes/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa
GTF=${LABSHR}/References/Annotations/hsapiens/gencode.v39.primary_assembly.annotation.gtf
OH=100

OUT=${SLURM_SUBMIT_DIR}/STAR

mkdir -p $OUT

$STAR --runMode genomeGenerate \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $OUT \
    --genomeFastaFiles $GENOME \
    --sjdbGTFfile $GTF \
    --sjdbOverhang $OH
