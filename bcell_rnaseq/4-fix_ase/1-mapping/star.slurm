#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=48gb
#SBATCH --time=24:00:00 
#SBATCH -p bch-compute 
#SBATCH --array=1-74
#SBATCH --job-name=STAR
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/STAR-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/STAR-%A-%a.err

source /programs/biogrids.shrc
STAR=/lab-share/IM-Gutierrez-e2/Public/vitor/software/STAR/STAR_2.7.11a/Linux_x86_64_static/STAR

# INPUT
META=${SLURM_SUBMIT_DIR}/metadata.tsv
DONORID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' $META )
SAMPLEID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' $META )
STIM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $3 }' $META )
FQ1=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $4 }' $META )
FQ2=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $5 }' $META )
FQDIR=${TEMP_WORK}/fastq/highinput 
FQ1qc=${FQDIR}/$( echo $FQ1 | sed 's|.fastq.gz|_val_1.fq.gz|' )
FQ2qc=${FQDIR}/$( echo $FQ2 | sed 's|.fastq.gz|_val_2.fq.gz|' )

LABSHR=/lab-share/IM-Gutierrez-e2/Public
VCF=${LABSHR}/vitor/lupus/bcell_rnaseq/4-fix_ase/0-genotypes/data/${DONORID}.vcf.gz

# INDEX
INDEX=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/indices/STAR

# OUTPUT
OUT=${TEMP_WORK}/bam/fixase/${SAMPLEID}_${STIM} 

# Run mapping
$STAR --runMode alignReads \
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

