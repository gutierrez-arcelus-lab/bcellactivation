#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=32gb
#SBATCH --time=48:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-5
#SBATCH --job-name=ASE
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ASE.%A.%a.out

source /programs/biogrids.shrc
export GATK_X=4.1.4.1
export PICARD_X=2.18.14
export SAMTOOLS_X=1.13
export QTLTOOLS_X=1.3.1

LABSHR=/lab-share/IM-Gutierrez-e2/Public
DATA=${LABSHR}/vitor/lupus/data 
RESULTS=${LABSHR}/vitor/lupus/bcell_bulk/results 

SAMPLELIST=${LABSHR}/vitor/lupus/bcell_bulk/bcell_samples.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )
ID=MG8989

BAMPREFIX=${TEMP_WORK}/star/pass_wasp/${SAMPLE} 
BAM=${BAMPREFIX}_Aligned.sortedByCoord.out.bam
GENOME=${DATA}/GRCh38.primary_assembly.genome.fa
VCF=${RESULTS}/gatk/${ID}.knownSNVs.biallelic.het.vcf
OUTDIR=${RESULTS}/ase
PREFIX=${OUTDIR}/${SAMPLE} 

mkdir -p $OUTDIR

# Extract uniquely mapped reads
samtools view -b -h -q 255 $BAM > ${BAMPREFIX}.uniq.bam

# Mark duplicates
picard AddOrReplaceReadGroups \
    I=${BAMPREFIX}.uniq.bam \
    O=${BAMPREFIX}.uniq.rg.bam \
    RGID=readGroupID${SLURM_ARRAY_TASK_ID} \
    RGLB=libraryID \
    RGPL=ILLUMINA \
    RGPU=HiSeq \
    RGSM=$ID \
    TMP_DIR=${TEMP_WORK}/temp

gatk MarkDuplicatesSpark \
    -I ${BAMPREFIX}.uniq.rg.bam \
    -O ${BAMPREFIX}.uniq.mkdups.bam \
    --tmp-dir ${TEMP_WORK}/temp \
    --spark-master local[$SLURM_CPUS_PER_TASK]

# Extract reads that passed WASP filter
samtools view -h ${BAMPREFIX}.uniq.mkdups.bam |\
    awk '$0 ~ /^@/ || $0 ~ /vW:i:1/' |\
    samtools view -Sb - > ${BAMPREFIX}.uniq.mkdups.wasp.bam

# Compute ASE
gatk ASEReadCounter \
    -R $GENOME \
    -I ${BAMPREFIX}.uniq.mkdups.wasp.bam \
    -V $VCF \
    -O ${PREFIX}.asereadcounter.txt \
    -min-depth 10 \
    --max-depth-per-sample 100

# QTLtools ase
GTF=${DATA}/gencode.v38.primary_assembly.annotation.gtf
GENOME2=${DATA}/GRCh38.primary_assembly.genome.fixnames.fa 

samtools sort -@ $SLURM_ARRAY_TASK_ID ${BAMPREFIX}.uniq.mkdups.wasp.bam > ${BAMPREFIX}.wasp.sorted.bam
samtools index -@ $SLURM_ARRAY_TASK_ID ${BAMPREFIX}.wasp.sorted.bam  

QTLtools ase \
    --bam ${BAMPREFIX}.wasp.sorted.bam \
    --vcf $VCF \
    --ind $ID \
    --mapq 255 \
    --fasta $GENOME2 \
    --gtf $GTF \
    --filter-duplicates \
    --cov 10 \
    --max-depth 100 \
    --out ${PREFIX}.qtltools

