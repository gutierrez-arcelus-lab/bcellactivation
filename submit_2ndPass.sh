#!/usr/benv bash

#PBS -l nodes=1:ppn=12
#PBS -l mem=48gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-2
#PBS -N STAR-2ndPass
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

SAMPLE=$( awk -v ARRID="$PBS_ARRAYID" 'FNR==ARRID { print $1 }' ./samples.txt )
FQDIR=/media/storage/genevol/geuvadis/fastq

FQ1=$FQDIR/${SAMPLE}_1.fastq.gz
FQ2=$FQDIR/${SAMPLE}_2.fastq.gz

INDEX=./star_index
SJDB=$( find ./pass_1st -name "*SJ*" | grep -f ./samples.txt)
OUT=./pass_2nd/${SAMPLE}

STAR --runMode alignReads \
    --runThreadN $PBS_NUM_PPN \
    --genomeDir $INDEX \
    --readFilesIn $FQ1 $FQ2 \
    --readFilesCommand zcat  \
    --outFilterMismatchNoverReadLmax 0.04 \
    --sjdbFileChrStartEnd $SJDB \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix ${OUT}_

