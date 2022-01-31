#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=48gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-5
#SBATCH --job-name=STARWASP
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/STAR-WASP.o%A_%a

source /programs/biogrids.shrc
export STAR_X=2.7.9a

LABSHR=/lab-share/IM-Gutierrez-e2/Public
SAMPLELIST=${LABSHR}/vitor/ase/bcell_bulk/bcell_samples.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )
ID=MG8989
FQDIR=${LABSHR}/B_cells/bulkTCpilot_1/34.198.31.178/210618_${ID}_fastq

FQ1=${FQDIR}/${SAMPLE}_R1_001.fastq.gz
FQ2=${FQDIR}/${SAMPLE}_R2_001.fastq.gz

INDEX=${LABSHR}/vitor/ase/indices/starindex
STAR2ND=${TEMP_WORK}/star/pass_2nd
SJDB=$( find $STAR2ND -name "*SJ*" | grep -f $SAMPLELIST)
VCF=${LABSHR}/vitor/ase/bcell_bulk/results/gatk/${ID}.knownSNVs.biallelic.het.vcf
OUTDIR=${LABSHR}/vitor/ase/bcell_bulk/results/star/pass_wasp
OUT=${OUTDIR}/${SAMPLE}

mkdir -p $OUTDIR

STAR --runMode alignReads \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $INDEX \
    --readFilesIn $FQ1 $FQ2 \
    --readFilesCommand zcat  \
    --outFilterMismatchNoverReadLmax 0.04 \
    --sjdbFileChrStartEnd $SJDB \
    --varVCFfile $VCF \
    --waspOutputMode SAMtag \
    --outSAMattributes vA vG vW \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${OUT}_

