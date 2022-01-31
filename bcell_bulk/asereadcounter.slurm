#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=48gb
#SBATCH --time=48:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-5
#SBATCH --job-name=ASE
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ASE.%A.%a.out

source /programs/biogrids.shrc
export SAMTOOLS_X=1.13
export PICARD_X=2.18.14
export GATK_X=4.1.4.1

LABSHR=/lab-share/IM-Gutierrez-e2/Public/vitor
SAMPLELIST=${LABSHR}/ase/bcell_bulk/bcell_samples.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )
ID=MG8989

BAMPREFIX=${LABSHR}/ase/bcell_bulk/results/star/pass_wasp/${SAMPLE} 
BAM=${BAMPREFIX}_Aligned.out.bam
GENOME=${LABSHR}/ase/data/GRCh38.primary_assembly.genome.fa
VCF=${LABSHR}/ase/bcell_bulk/results/gatk/${ID}.knownSNVs.biallelic.het.vcf
OUTDIR=${LABSHR}/ase/bcell_bulk/results/gatk
PREFIX=${OUTDIR}/${SAMPLE} 
BAMTMP=${TEMP_WORK}/gatk/${SAMPLE}


# Extract uniquely mapped reads
samtools view -b -h -q 255 $BAM > ${BAMTMP}.uniq.bam

# Mark duplicates
picard AddOrReplaceReadGroups \
    I=${BAMTMP}.uniq.bam \
    O=${BAMTMP}.uniq.rg.bam \
    RGID=readGroupID${SLURM_ARRAY_TASK_ID} \
    RGLB=libraryID \
    RGPL=ILLUMINA \
    RGPU=HiSeq \
    RGSM=$SAMPLE

gatk MarkDuplicatesSpark \
    -I ${BAMTMP}.uniq.rg.bam \
    -O ${BAMTMP}.uniq.mkdups.bam \
    --spark-master local[$SLURM_CPUS_PER_TASK]

# Extract reads that passed WASP filter
samtools view -h ${BAMTMP}.uniq.mkdups.bam |\
    awk '$0 ~ /^@/ || $0 ~ /vW:i:1/' |\
    samtools view -Sb - > ${BAMTMP}.uniq.mkdups.wasp.bam

# Compute ASE
if [[ ! -f "$VCF".idx  ]]; then
    gatk IndexFeatureFile -I $VCF
fi

gatk ASEReadCounter \
    -R $GENOME \
    -I ${BAMTPM}.uniq.mkdups.wasp.bam \
    -V $VCF \
    -O ${PREFIX}.ase.counts \
    -min-depth 8
