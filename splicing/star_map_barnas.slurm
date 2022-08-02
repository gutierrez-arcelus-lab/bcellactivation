#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=48gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute-pe
#SBATCH --array=1-13
#SBATCH --job-name=STAR
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/STAR-Barnas.o%A.%a

source /programs/biogrids.shrc
export STAR_X=2.7.9a
export REGTOOLS_X=0.5.1
export SAMTOOLS_X=1.13

LABSHR=/lab-share/IM-Gutierrez-e2/Public 
SAMPLELIST=${LABSHR}/External_datasets/Barnas/SRR_Acc_List_DN.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )

FQDIR=${LABSHR}/External_datasets/Barnas/fastq
FQ=${FQDIR}/${SAMPLE}.fastq.gz

INDEX=${SLURM_SUBMIT_DIR}/starindex_101bp
OUTDIR=${SLURM_SUBMIT_DIR}/mapping/barnas
OUT=${OUTDIR}/${SAMPLE}

mkdir -p $OUTDIR

# Map reads
STAR --runMode alignReads \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $INDEX \
    --readFilesIn $FQ \
    --readFilesCommand zcat  \
    --outFilterMismatchNoverReadLmax 0.04 \
    --twopassMode Basic \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${OUT}_

rm -r ${OUT}__*

# Extract junctions
BAM=${OUT}_Aligned.sortedByCoord.out.bam
samtools index -@ $SLURM_CPUS_PER_TASK $BAM
regtools junctions extract -s 0 -a 8 -m 50 -M 500000 $BAM -o ${OUT}.junc

