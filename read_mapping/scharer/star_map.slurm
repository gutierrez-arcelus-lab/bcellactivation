#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-81
#SBATCH --job-name=STAR
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/STAR-Scharer.o%A.%a

source /programs/biogrids.shrc
export STAR_X=2.7.9a
export REGTOOLS_X=0.5.1
export SAMTOOLS_X=1.13

LABSHR=/lab-share/IM-Gutierrez-e2/Public 
SAMPLELIST=${LABSHR}/External_datasets/Scharer/RNAseq/metadata/SRR_Acc_List.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )

FQDIR=${LABSHR}/External_datasets/Scharer/RNAseq/fastq
FQ1=$FQDIR/${SAMPLE}_1.fastq.gz
FQ2=$FQDIR/${SAMPLE}_2.fastq.gz

INDEX=${LABSHR}/vitor/lupus/indices/starindex_50
OUTDIR=${SLURM_SUBMIT_DIR}/mapping
OUT=${OUTDIR}/${SAMPLE}

mkdir -p $OUTDIR

# Map reads
STAR --runMode alignReads \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $INDEX \
    --readFilesIn $FQ1 $FQ2 \
    --readFilesCommand zcat  \
    --outFilterMismatchNoverReadLmax 0.04 \
    --twopassMode Basic \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${OUT}_

# Extract junctions
BAM=${OUT}_Aligned.sortedByCoord.out.bam
samtools index -@ $SLURM_CPUS_PER_TASK $BAM
samtools view -@ $SLURM_CPUS_PER_TASK -b -h -q 255 $BAM `(seq 1 22; echo X) | sed 's/^/chr/'` > ${OUT}.uniq.bam
samtools index -@ $SLURM_CPUS_PER_TASK ${OUT}.uniq.bam
regtools junctions extract -s 0 -a 8 -m 50 -M 500000 ${OUT}.uniq.bam -o ${OUT}.junc

rm ${OUT}.uniq.bam*
rm -r ${OUT}__*
