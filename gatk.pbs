#!/usr/bin/env bash

#PBS -l nodes=1:ppn=16
#PBS -l mem=64gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-2
#PBS -N GATK
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

SAMPLE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $1 }' ./samples.txt )

BAM=./output/pass_2nd/${SAMPLE}_Aligned.out.bam
OUT=./output/gatk/${SAMPLE}_rg_added_sorted.bam 

PICARD=/home/vitor/Libraries/picard.jar
GATK=/home/vitor/Libraries/gatk-4.2.2.0/gatk

# Add read groups, sort, mark duplicates, and create index

java -Xmx8g -jar $PICARD AddOrReplaceReadGroups \
    -I $BAM \
    -O $OUT \
    -SO coordinate \
    -RGID readGroupID \
    -RGLB libraryID \
    -RGPL illumina \
    -RGPU HiSeq \
    -RGSM sampleID 

OUT2=./output/gatk/${SAMPLE}_dedupped.bam

java -Xmx8g -jar $PICARD MarkDuplicates \
    -I $OUT \
    -O $OUT2 \
    -CREATE_INDEX true \
    -VALIDATION_STRINGENCY SILENT \
    -M ./output/gatk/${SAMPLE}_output_metrics

# Split'N'Trim and reassign mapping qualities

GENOME=/home/vitor/gencode/GRCh38.primary_assembly.genome.fa
OUT3=./output/gatk/${SAMPLE}_split.bam

if [[ ! -f "$GENOME".fai ]]; then
    samtools faidx $GENOME
fi

if [[ ! -f "$GENOME".dict ]]; then
    $GATK CreateSequenceDictionary -R $GENOME
fi

# The following options cannot be used with latest GATK
# However, 255->60 quality transformer seems to be on by default
#
#-RF ReassignOneMappingQuality \
#-RMQF 255 \
#-RMQT 60 \
#-U ALLOW_N_CIGAR_READS

$GATK --java-options "-Xmx8g" SplitNCigarReads \
    -R $GENOME \
    -I $OUT2 \
    -O $OUT3

# Base Recalibration

# I could not run this step because
# need to decide what to use as "knownSites"

$GATK --java-options "-Xmx8g" BaseRecalibrator \
    -I $OUT3 \
    -R $GENOME \
    -knownSites /humgen/gsa-hpprojects/GATK/bundle/current/hg19/dbsnp_138.hg19.vcf \
    -knownSites /humgen/gsa-hpprojects/GATK/bundle/current/hg19/1000G_phase1.snps.high_confidence.hg19.vcf \
    -knownSites /humgen/gsa-hpprojects/GATK/bundle/current/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf \
    -o recal_data.table \
    -nct 4

$GATK --java-options "-Xmx64g" PrintReads \
    -I $OUT3 \
    -R $GENOME \
    -BQSR recal_data.table \
    -o split-BQSR.bam \
    -nct 4
