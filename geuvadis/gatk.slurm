#! /usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1
#SBATCH --job-name=GATK
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/GATK.%A.%a.out

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12
export PICARD_X=2.18.14


SAMPLELIST=/home/ch229163/ase/samples.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )

LABSHARE=/lab-share/IM-Gutierrez-e2/Public/vitor
BAM=${LABSHARE}/ase/geuvadis/results/star/pass_2nd/${SAMPLE}_Aligned.out.bam
OUTDIR=${LABSHARE}/geuvadis/results/gatk
PREFIX=${OUTDIR}/${SAMPLE} 
OUT=${PREFIX}_rg.bam 

mkdir -p $OUTDIR

# Add a read group to the BAM file
picard AddOrReplaceReadGroups \
    -I $BAM \
    -O $OUT \
    -RGID readGroupID \
    -RGLB libraryID \
    -RGPL ILLUMINA \
    -RGPU HiSeq \
    -RGSM $SAMPLE 

# Mark duplicated reads
OUT2=${PREFIX}_deduplicated.bam 

gatk MarkDuplicatesSpark \
    -I $OUT \
    -O $OUT2 \
    -M ${PREFIX}_metrics.txt \
    --spark-master local[$SLURM_CPUS_PER_TASK]


# Split'N'Trim and reassign mapping qualities
GENOME=$HOME/ase/data/GRCh38.primary_assembly.genome.fa
OUT3=${PREFIX}_splitn.bam

if [[! -f "$GENOME" ]]; then
    zcat ${GENOME}.gz > $GENOME
fi

if [[ ! -f "$GENOME".fai ]]; then
    samtools faidx $GENOME
fi

if [[ ! -f "$GENOME".dict ]]; then
    $GATK CreateSequenceDictionary -R $GENOME
fi

gatk SplitNCigarReads \
    -R $GENOME \
    -I $OUT2 \
    -O $OUT3

# Base Recalibration
gatk BaseRecalibrator \
    -R $GENOME \
    -I $OUT3 \
    -O ${PREFIX}_bqsr.report \
    --known-sites $HOME/ase/data/dbsnp_155.hg38.vcf.gz

# Apply recallibration
OUT4=${PREFIX}_recal.bam 

gatk ApplyBQSR \
    -R $GENOME \
    -I $OUT3 \
    -O $OUT4 \
    --bqsr-recal-file ${PREFIX}_bqsr.report  

# HaplotypeCaller
VCF=${PREFIX}.vcf.gz

gatk HaplotypeCaller \
    -R $GENOME \
    -I $OUT4 \
    -O $VCF \
    --dont-use-soft-clipped-bases true \
    -stand-call-conf 10 \	
    --native-pair-hmm-threads $SLURM_CPUS_PER_TASK

FILTVCF=${PREFIX}.filtered.vcf.gz
LOG=${PREFIX}.snpCluster.log 

gatk VariantFiltration \
    -R $GENOME \
    -V $VCF \
    -O $FILTVCF \
    -window 35 -cluster 3 &> $LOG 

bcftools norm -m -both $FILTVCF -O z -o ${PREFIX}.filtered.denorm.vcf.gz
