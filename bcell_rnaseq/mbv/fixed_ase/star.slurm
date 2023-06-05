#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=48gb
#SBATCH --time=24:00:00 
#SBATCH -p bch-compute 
#SBATCH --array=1-17
#SBATCH --job-name=STAR
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/STAR-%A-%a

source /programs/biogrids.shrc
export STAR_X=2.7.9a

# INPUT
FQDIR=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/QC 
META=${SLURM_SUBMIT_DIR}/data/fix_metadata.tsv
SUBJECTID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' $META )
SAMPLEID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' $META )
STIM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $3 }' $META )
FQ1=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $4 }' $META )
FQ2=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $5 }' $META )
FQ1qc=${FQDIR}/$( echo $FQ1 | sed 's|.fastq.gz|_val_1.fq.gz|' )
FQ2qc=${FQDIR}/$( echo $FQ2 | sed 's|.fastq.gz|_val_2.fq.gz|' )

VCF=${SLURM_SUBMIT_DIR}/data/${SUBJECTID}.vcf.gz

# INDEX
INDEX=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/indices/starindex_100

# OUTPUT
OUT=${SLURM_SUBMIT_DIR}/mapping/${SAMPLEID}_${STIM}

# EXECUTE
STAR --runMode alignReads \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $INDEX \
    --readFilesIn $FQ1qc $FQ2qc \
    --readFilesCommand zcat  \
    --outFilterMismatchNoverReadLmax 0.04 \
    --twopassMode Basic \
    --varVCFfile <(zcat $VCF) \
    --waspOutputMode SAMtag \
    --outSAMattributes vA vG vW \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${OUT}_

