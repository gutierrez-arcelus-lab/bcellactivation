#! /usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-23
#SBATCH --job-name=HapCall
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/HaplotypeCaller.%A.%a.out

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12
export GATK_X=4.1.4.1
export SAMTOOLS_X=1.13

# Chromosome name
if [ "$SLURM_ARRAY_TASK_ID" -lt 23 ]; then
    
    CHR=chr${SLURM_ARRAY_TASK_ID}

elif [ "$SLURM_ARRAY_TASK_ID" -eq 23 ]; then
    
    CHR=chrX

else
    echo "ERROR in chr spec"
fi

# IO
LABSHR=/lab-share/IM-Gutierrez-e2/Public
ID=MG8989
SAMPLES=( $(cat ${LABSHR}/vitor/lupus/bcell_bulk/bcell_samples.txt) )
BAMS=( "${SAMPLES[@]/#/$TEMP_WORK/gatk/}" )
BAMS=( "${BAMS[@]/%/_recal.bam}" )
GENOME=${LABSHR}/vitor/lupus/data/GRCh38.primary_assembly.genome.fa
OUTDIR=${TEMP_WORK}/gatk
PREFIX=${OUTDIR}/${ID}.${CHR}
MERGED=${PREFIX}.merged.bam
VCF=${PREFIX}.vcf.gz

# Merge stim conditions into single BAM
samtools merge -@ $SLURM_CPUS_PER_TASK -R $CHR -o $MERGED "${BAMS[@]}"
samtools index -@ $SLURM_CPUS_PER_TASK $MERGED

# HaplotypeCaller
gatk HaplotypeCaller \
    -R $GENOME \
    -I $MERGED \
    -O $VCF \
    --dont-use-soft-clipped-bases true \
    -stand-call-conf 10 \
    --native-pair-hmm-threads $SLURM_CPUS_PER_TASK \
    --tmp-dir ${TEMP_WORK}/temp

# Variant filtering
FILTVCF=${PREFIX}.filtered.vcf.gz
LOG=${PREFIX}.snpCluster.log 

gatk VariantFiltration \
    -R $GENOME \
    -V $VCF \
    -O $FILTVCF \
    --tmp-dir ${TEMP_WORK}/temp \
    -window 35 -cluster 3 &> $LOG 

# Save only SNVs
FILTVCF_SNV=${PREFIX}.filtered.biallelic.het.SNVs.vcf.gz 

bcftools view -m2 -M2 -v snps -O z -o $FILTVCF_SNV $FILTVCF 

# Extract only SNVs reported by dbSNP
DBSNP=${LABSHR}/vitor/lupus/data/dbsnp_155.hg38pri.vcf.gz
POS=${PREFIX}.pos
DBSNP_VCF=${PREFIX}.dbSNP.vcf.gz

# Extract dbSNP SNV variants for the current sample
zcat $FILTVCF_SNV |\
    grep -v "^#" |\
    awk '{ print $1"\t"$2 }' |\
    uniq > $POS

bcftools view -R $POS -v snps -O z -o $DBSNP_VCF $DBSNP

# Now refilter the sample VCF for the SNV variants present in dbSNP
FINAL_POS=${PREFIX}.final.pos

zcat $DBSNP_VCF |\
    grep -v "^#" |\
    awk '{ print $1"\t"$2 }' |\
    uniq > $FINAL_POS

FINAL_VCF=${PREFIX}.knownSNVs.biallelic.het.vcf

tabix -p vcf $FILTVCF_SNV  
bcftools view -R $FINAL_POS -O v -o $FINAL_VCF $FILTVCF_SNV 

