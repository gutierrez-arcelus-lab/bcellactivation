#! /usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-5
#SBATCH --job-name=GATK
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/GATK.%A.%a.out

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12
export PICARD_X=2.18.14
export GATK_X=4.1.4.1
export SAMTOOLS_X=1.13

LABSHR=/lab-share/IM-Gutierrez-e2/Public
SAMPLELIST=${LABSHR}/vitor/lupus/bcell_bulk/bcell_samples.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )
ID=MG8989
BAM=${TEMP_WORK}/star/pass_2nd/${SAMPLE}_Aligned.out.bam
BAMUNIQ=${TEMP_WORK}/star/pass_2nd/${SAMPLE}_uniq.bam
OUTDIR=${TEMP_WORK}/gatk
PREFIX=${OUTDIR}/${SAMPLE}

mkdir -p $OUTDIR

samtools view -b -h -q 255 $BAM > $BAMUNIQ

# Add a read group to the BAM file
OUT=${PREFIX}_rg.bam 

picard AddOrReplaceReadGroups \
    I=$BAMUNIQ \
    O=$OUT \
    RGID=readGroup${SLURM_ARRAY_TASK_ID} \
    RGLB=libraryID \
    RGPL=ILLUMINA \
    RGPU=HiSeq \
    RGSM=$ID \
    --TMP_DIR=$TEMP_WORK/temp \

# Mark duplicated reads
OUT2=${PREFIX}_mkdups.bam 

gatk MarkDuplicatesSpark \
    -I $OUT \
    -O $OUT2 \
    -M ${PREFIX}_metrics.txt \
    --tmp-dir $TEMP_WORK/temp \
    --spark-master local[$SLURM_CPUS_PER_TASK]


# Split'N'Trim and reassign mapping qualities
GENOME=${LABSHR}/vitor/lupus/data/GRCh38.primary_assembly.genome.fa
OUT3=${PREFIX}_splitn.bam

if [[ ! -f "$GENOME" ]]; then
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
    -O $OUT3 \
    --tmp-dir $TEMP_WORK/temp

# Base Recalibration
DBSNP=${LABSHR}/vitor/lupus/data/dbsnp_155.hg38pri.vcf.gz

gatk BaseRecalibrator \
    -R $GENOME \
    -I $OUT3 \
    -O ${PREFIX}_bqsr.report \
    --known-sites $DBSNP \
    --tmp-dir $TEMP_WORK/temp

# Apply recallibration
OUT4=${PREFIX}_recal.bam 

gatk ApplyBQSR \
    -R $GENOME \
    -I $OUT3 \
    -O $OUT4 \
    --bqsr-recal-file ${PREFIX}_bqsr.report \
    --tmp-dir $TEMP_WORK/temp

