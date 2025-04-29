#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=48gb
#SBATCH --time=48:00:00 
#SBATCH -p bch-compute 
#SBATCH --array=1-74
#SBATCH --job-name=STAR
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/STAR-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/STAR-%A-%a.err

source /programs/biogrids.shrc
export STAR_X=2.7.11a
export REGTOOLS_X=0.5.1
export SAMTOOLS_X=1.13

# INPUT
META=${SLURM_SUBMIT_DIR}/metadata_qced.tsv
SAMPLEID=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' )
FQ1=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' ) 
FQ2=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $3 }' ) 

# INDEX
INDEX=${SLURM_SUBMIT_DIR}/index/STAR

# OUTPUT
OUT=${SLURM_SUBMIT_DIR}/bam/${SAMPLEID} 

# Run mapping
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