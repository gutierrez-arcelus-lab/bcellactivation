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
#SBATCH -o /temp_work/ch229163/log/ASE-%A-%a

source /programs/biogrids.shrc
export GATK_X=4.1.4.1
export PICARD_X=2.18.14
export SAMTOOLS_X=1.13

CPUS=$SLURM_CPUS_PER_TASK

# INPUT
GENOME=/lab-share/IM-Gutierrez-e2/Public/References/Genomes/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa

META=${SLURM_SUBMIT_DIR}/arrayspec_ase.tsv
SUBJECTID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $META )
SAMPLEID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $META )
STIM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $3 }' $META )
MGBID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $4 }' $META )

BAM=${SLURM_SUBMIT_DIR}/mapping/${SAMPLEID}_${STIM}_Aligned.sortedByCoord.out.bam
VCF=${SLURM_SUBMIT_DIR}/data/${SUBJECTID}.vcf.gz

PREFIX=${SLURM_SUBMIT_DIR}/ase/${SAMPLEID}_${STIM}
TEMP=${TEMP_WORK}/${SAMPLEID}_${STIM} 

# Extract uniquely mapped reads
samtools view --threads $CPUS -b -h -q 255 -o ${TEMP}.uniq.bam $BAM 

# Mark duplicates
picard AddOrReplaceReadGroups \
    I=${TEMP}.uniq.bam \
    O=${TEMP}.uniq.rg.bam \
    RGID=readGroupID${SLURM_ARRAY_TASK_ID} \
    RGLB=libraryID \
    RGPL=ILLUMINA \
    RGPU=NovaSeq \
    RGSM=$MGBID \
    TMP_DIR=${TEMP_WORK}/temp

gatk MarkDuplicatesSpark \
    -I ${TEMP}.uniq.rg.bam \
    -O ${TEMP}.uniq.mkdups.bam \
    --tmp-dir ${TEMP_WORK}/temp \
    --spark-master local[$CPUS]

# Extract reads that passed WASP filter
samtools view -h ${TEMP}.uniq.mkdups.bam |\
    awk '$0 ~ /^@/ || $0 ~ /vW:i:1/' |\
    samtools view -Sb - > ${TEMP}.uniq.mkdups.wasp.bam

# Compute ASE
gatk ASEReadCounter \
    -R $GENOME \
    -I ${TEMP}.uniq.mkdups.wasp.bam \
    -V $VCF \
    -O ${PREFIX}.asereadcounter.txt \
    -min-depth 10 \
    --max-depth-per-sample 100
