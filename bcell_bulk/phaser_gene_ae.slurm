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
SAMPLELIST=( $( cat ${LABSHR}/vitor/lupus/bcell_bulk/bcell_samples.txt ) )
RESULTS=${LABSHR}/vitor/lupus/bcell_bulk/results 

ID=MG8989
BAMS=( "${SAMPLELIST[@]/#/$TEMP_WORK/star/pass_wasp/}" )
BAMS=( "${BAMS[@]/%/.wasp.sorted.bam}" )
VCF=${RESULTS}/phaser/MG8989_merge.vcf.gz
OUTDIR=${RESULTS}/phaser
PREFIX=${OUTDIR}/${ID}

python $HOME/phaser/phaser/phaser.py \
    --vcf $VCF \
    --bam ${BAMS[0]},${BAMS[1]},${BAMS[2]},${BAMS[3]},${BAMS[4]} \
    --sample $ID \
    --paired_end 1 --mapq 255 --baseq 10 --pass_only 0 \
    --threads $SLURM_CPUS_PER_TASK \
    --o $PREFIX

python $HOME/phaser/phaser_gene_ae/phaser_gene_ae.py \
    --haplotypic_counts ${PREFIX}.haplotypic_counts.txt \
    --features $HOME/lupus/data/gencodev38_genes.bed \
    --o ${PREFIX}_gene_ae.txt
