#! /usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1
#SBATCH --job-name=GATK
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/GATK.%A.%a.out

source /programs/biogrids.shrc
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


## Split'N'Trim and reassign mapping qualities
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
    --known-sites $HOME/ase/data/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf \
    --known-sites $HOME/ase/data/GCF_000001405.39.gz \
    --known-sites $HOME/ase/data/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites $HOME/ase/data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 

## Next: ApplyBQSR replaces PrintReads in GATK4
#
#$GATK --java-options "-Xmx64g" ApplyBQSR \
#    -R $GENOME \
#    -I $OUT3 \
#    --bqsr-recal-file ./output/gatk/${SAMPLE}_recal_data.table \
#    -O ./output/gatk/${SAMPLE}_split_BQSR.bam
#
