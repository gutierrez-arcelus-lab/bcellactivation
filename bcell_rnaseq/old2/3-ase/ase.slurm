#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=48gb
#SBATCH --time=24:00:00 
#SBATCH -p bch-compute 
#SBATCH --array=2-75
#SBATCH --job-name=ASE
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ASE-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/ASE-%A-%a.err

source /programs/biogrids.shrc
export STAR_X=2.7.11a
export GATK_X=4.1.4.1
export PICARD_X=2.18.14
export SAMTOOLS_X=1.13

LABSHR=/lab-share/IM-Gutierrez-e2/Public

# INPUT
META=${LABSHR}/vitor/lupus/bcell_rnaseq/2-mbv/matching_results.tsv
RNASAMPLEID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' $META )
STIM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' $META )
VCFDONORID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $3 }' $META )
RNAREP=$(echo $RNASAMPLEID | cut -d'_' -f2)

FQ1=${TEMP_WORK}/fastq/highinput_merged/${RNASAMPLEID}_${STIM}_R1_val_1.fq.gz
FQ2=${TEMP_WORK}/fastq/highinput_merged/${RNASAMPLEID}_${STIM}_R2_val_2.fq.gz

VCF=${SLURM_SUBMIT_DIR}/0-genotypes/data/${VCFDONORID}.vcf.gz

# INDEX
INDEX=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/indices/STAR

# OUTPUT
OUTDIR=${TEMP_WORK}/bam/ase_v2/${VCFDONORID}_${RNAREP}_${STIM} 
mkdir -p $OUTDIR

# Run read mapping
STAR --runMode alignReads \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $INDEX \
    --readFilesIn $FQ1 $FQ2 \
    --readFilesCommand zcat  \
    --outFilterMismatchNoverReadLmax 0.04 \
    --twopassMode Basic \
    --varVCFfile <(zcat $VCF) \
    --waspOutputMode SAMtag \
    --outSAMattributes vA vG vW \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${OUTDIR}/

# Extract uniquely mapped reads
samtools view --threads $SLURM_CPUS_PER_TASK -b -h -q 255 \
    -o ${OUTDIR}/uniq.bam \
    ${OUTDIR}/Aligned.sortedByCoord.out.bam

# Mark duplicates
picard AddOrReplaceReadGroups \
    I=${OUTDIR}/uniq.bam \
    O=${OUTDIR}/uniq_rg.bam \
    RGID=RG${SLURM_ARRAY_TASK_ID} \
    RGLB=libraryID \
    RGPL=ILLUMINA \
    RGPU=NovaSeq \
    RGSM=${VCFDONORID} \
    TMP_DIR=${OUTDIR}

gatk MarkDuplicatesSpark \
    -I ${OUTDIR}/uniq_rg.bam \
    -O ${OUTDIR}/uniq_mkdups.bam \
    --tmp-dir $OUTDIR \
    --spark-master local[$SLURM_CPUS_PER_TASK]

# Extract reads that pass WASP
samtools view -h ${OUTDIR}/uniq_mkdups.bam |\
    awk '$0 ~ /^@/ || $0 ~ /vW:i:1/' |\
    samtools view -Sb - > ${OUTDIR}/uniq_mkdups_wasp.bam

# ASE
GENOME=/lab-share/IM-Gutierrez-e2/Public/References/Genomes/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa

gatk ASEReadCounter \
    -R $GENOME \
    -I ${OUTDIR}/uniq_mkdups_wasp.bam \
    -V $VCF \
    -O ${SLURM_SUBMIT_DIR}/results/${VCFDONORID}_${RNAREP}_${STIM}_ase.txt 
