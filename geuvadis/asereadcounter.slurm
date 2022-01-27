#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1
#SBATCH --job-name=ASE
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ASE.%A.%a.out

source /programs/biogrids.shrc
export SAMTOOLS_X=1.14
export PICARD_X=2.18.14
export GATK_X=4.1.4.1

LABSHARE=/lab-share/IM-Gutierrez-e2/Public/vitor
SAMPLELIST=${LABSHARE}/ase/geuvadis/samples.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )

BAMPREFIX=${LABSHARE}/ase/geuvadis/results/star/pass_wasp/${SAMPLE} 
BAM=${BAMPREFIX}_Aligned.out.bam
GENOME=${LABSHARE}/ase/data/GRCh38.primary_assembly.genome.fa
VCF=${LABSHARE}/ase/geuvadis/results/gatk/${SAMPLE}.knownSNVs.biallelic.het.vcf
OUTDIR=${LABSHARE}/ase/geuvadis/results/gatk
PREFIX=${OUTDIR}/${SAMPLE} 


# Extract uniquely mapped reads
samtools view -b -h -q 255 $BAM > ${BAMPREFIX}.uniq.bam

# Mark duplicates
picard AddOrReplaceReadGroups \
    I=${BAMPREFIX}.uniq.bam \
    O=${BAMPREFIX}.uniq.rg.bam \
    RGID=readGroupID \
    RGLB=libraryID \
    RGPL=ILLUMINA \
    RGPU=HiSeq \
    RGSM=$SAMPLE

gatk MarkDuplicatesSpark \
    -I ${BAMPREFIX}.uniq.rg.bam \
    -O ${BAMPREFIX}.uniq.mkdups.bam \
    --spark-master local[$SLURM_CPUS_PER_TASK]

# Extract reads that passed WASP filter
samtools view -h ${BAMPREFIX}.uniq.mkdups.bam |\
    awk '$0 ~ /^@/ || $0 ~ /vW:i:1/' |\
    samtools view -Sb - > ${BAMPREFIX}.uniq.mkdups.wasp.bam

# Compute ASE
if [[ ! -f "$VCF".idx  ]]; then
    gatk IndexFeatureFile -I $VCF
fi

gatk ASEReadCounter \
    -R $GENOME \
    -I ${BAMPREFIX}.uniq.mkdups.bam \
    -V $VCF \
    -O ${PREFIX}.ase.counts \
    -min-depth 8

gatk ASEReadCounter \
    -R $GENOME \
    -I ${BAMPREFIX}.uniq.mkdups.wasp.bam \
    -V $VCF \
    -O ${PREFIX}.ase.wasp.counts \
    -min-depth 8
