#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=36gb
#SBATCH -p bch-compute
#SBATCH --job-name=STARidx
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBACTH --mail-type=END,FAIL
#SBATCH -o logs/STAR-index.o%A-%a.out
#SBATCH -e logs/STAR-index.o%A-%a.err

source /programs/biogrids.shrc
export STAR_X=2.7.11a

cd "$SLURM_SUBMIT_DIR"

GENOME="./data/GRCh38.primary_assembly.genome.fa"
GTF="./data/gencode.v38.primary_assembly.annotation.gtf"
OH=100

OUT="./index"

mkdir -p $OUT

STAR --runMode genomeGenerate \
    --runThreadN "$SLURM_CPUS_PER_TASK" \
    --genomeDir "$OUT" \
    --genomeFastaFiles "$GENOME" \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang "$OH"
