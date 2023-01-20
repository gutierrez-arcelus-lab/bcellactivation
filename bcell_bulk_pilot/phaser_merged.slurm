#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=PHASER
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/PHASER.%j.out

source /programs/biogrids.shrc

export SAMTOOLS_X=1.13
export BEDTOOLS_X=2.29.0

LABSHR=/lab-share/IM-Gutierrez-e2/Public
RESULTS=${LABSHR}/vitor/lupus/bcell_bulk/results 

SAMPLES=( $(cat ${LABSHR}/vitor/lupus/bcell_bulk/bcell_samples.txt) )
ID=MG8989
BAMS=( "${SAMPLES[@]/#/$TEMP_WORK/star/pass_wasp/}" )
BAMS=( "${BAMS[@]/%/.wasp.sorted.bam}" )
MERGE=${TEMP_WORK}/star/pass_wasp/${ID}.wasp.merged.bam

samtools merge -@ $SLURM_CPUS_PER_TASK -o $MERGE "${BAMS[@]}"
samtools index -@ $SLURM_CPUS_PER_TASK $MERGE

VCF=${RESULTS}/gatk/${ID}.knownSNVs.biallelic.het.vcf.gz
OUTDIR=${RESULTS}/phaser
PREFIX=${OUTDIR}/${ID}_merge 

mkdir -p $OUTDIR

python $HOME/phaser/phaser/phaser.py \
    --vcf $VCF \
    --bam $MERGE \
    --sample $ID \
    --paired_end 1 --mapq 255 --baseq 10 --pass_only 0 \
    --gw_phase_vcf 1 \
    --threads $SLURM_CPUS_PER_TASK \
    --o $PREFIX

