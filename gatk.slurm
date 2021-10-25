#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1
#SBATCH --job-name=GATK
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /home/ch229163/ase/log/GATK.%A.%a.out

source /programs/biogrids.shrc
export PICARD_X=2.18.14


SAMPLELIST=/home/ch229163/ase/samples.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )

LABSHARE=/lab-share/IM-Gutierrez-e2/Public/vitor
BAM=$LABSHARE/pass_2nd/${SAMPLE}_Aligned.out.bam
OUTDIR=$LABSHARE/gatk
OUT=$OUTDIR/${SAMPLE}_rg_added_sorted.bam 

PICARD=/home/vitor/Libraries/picard.jar
GATK=/home/vitor/Libraries/gatk-4.2.2.0/gatk

# Add a single read group to all reads in a file
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


$GATK --java-options "-Xmx8g" SplitNCigarReads \
    -R $GENOME \
    -I $OUT2 \
    -O $OUT3

# The following options cannot be used with latest GATK
# However, 255->60 quality transformer seems to be on by default
#
#-RF ReassignOneMappingQuality \
#-RMQF 255 \
#-RMQT 60 \
#-U ALLOW_N_CIGAR_READS

# Base Recalibration

# To check:
# the knownSites files downloaded from GATK
# have contig IDs that are different from the index

$GATK --java-options "-Xmx8g" BaseRecalibrator \
    -R $GENOME \
    -I $OUT3 \
    --known-sites ./data/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites ./data/1000G_phase1.snps.high_confidence.hg38.vcf \
    --known-sites ./data/Mills_and_1000G_gold_standard.indels.hg38.vcf \
    -O ./output/gatk/${SAMPLE}_recal_data.table

# There is no -nct option in GATK4
#-nct 4
# Alternative would be to use SPARK,
# but I couldn't make it work

# Next: ApplyBQSR replaces PrintReads in GATK4

$GATK --java-options "-Xmx64g" ApplyBQSR \
    -R $GENOME \
    -I $OUT3 \
    --bqsr-recal-file ./output/gatk/${SAMPLE}_recal_data.table \
    -O ./output/gatk/${SAMPLE}_split_BQSR.bam

