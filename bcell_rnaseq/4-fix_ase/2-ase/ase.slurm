#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=36gb
#SBATCH --time=48:00:00 
#SBATCH -p bch-compute 
#SBATCH --array=1-74
#SBATCH --job-name=ASE
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ASE-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/ASE-%A-%a.err

source /programs/biogrids.shrc
export GATK_X=4.1.4.1
export PICARD_X=2.18.14
export SAMTOOLS_X=1.13

CPUS=$SLURM_CPUS_PER_TASK

# INPUT
## GENOME
LABSHR=/lab-share/IM-Gutierrez-e2/Public 
GENOME=${LABSHR}/References/Genomes/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa

## SAMPLE METADATA
SPEC=${SLURM_SUBMIT_DIR}/array_spec.tsv
DONORID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SPEC )
SAMPLEID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $SPEC )
STIM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $3 }' $SPEC )
MGBID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $4 }' $SPEC )

## BAM and VCF
PREFIX=${TEMP_WORK}/bam/fixase/${SAMPLEID}_${STIM} 
BAM=${PREFIX}_Aligned.sortedByCoord.out.bam
VCF=${LABSHR}/vitor/lupus/bcell_rnaseq/4-fix_ase/0-genotypes/data/${DONORID}.vcf.gz

# OUTPUT
OUT=${SLURM_SUBMIT_DIR}/results 
mkdir -p $OUT

# Create a temporary directory for each sample and delete it at the end
TEMP=${TEMP_WORK}/temp/gatk/${SAMPLEID}_${STIM} 
mkdir -p $TEMP

# Extract uniquely mapped reads
samtools view --threads $CPUS -b -h -q 255 -o ${TEMP}/uniq.bam $BAM 

# Mark duplicates
picard AddOrReplaceReadGroups \
    I=${TEMP}/uniq.bam \
    O=${TEMP}/uniq.rg.bam \
    RGID=readGroupID${SLURM_ARRAY_TASK_ID} \
    RGLB=libraryID \
    RGPL=ILLUMINA \
    RGPU=NovaSeq \
    RGSM=$MGBID \
    TMP_DIR=$TEMP

gatk MarkDuplicatesSpark \
    -I ${TEMP}/uniq.rg.bam \
    -O ${TEMP}/uniq.mkdups.bam \
    --tmp-dir $TEMP \
    --spark-master local[$CPUS]

# Extract reads that passed WASP filter.
samtools view -h ${TEMP}/uniq.mkdups.bam |\
    awk '$0 ~ /^@/ || $0 ~ /vW:i:1/' |\
    samtools view -Sb - > ${PREFIX}.uniq.mkdups.wasp.bam

# Compute ASE
# filters out duplicate reads by default
gatk ASEReadCounter \
    -R $GENOME \
    -I ${PREFIX}.uniq.mkdups.wasp.bam \
    -V $VCF \
    -O ${OUT}/${SAMPLEID}_${STIM}.asereadcounter.txt
    
#-min-depth 1
#--max-depth-per-sample 60

# Remove TEMP directory
#rm -r $TEMP
