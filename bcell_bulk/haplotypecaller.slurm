#! /usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64gb
#SBATCH --time=72:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=HapCall
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/HaplotypeCaller.%A.%a.out

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12
export GATK_X=4.1.4.1
export SAMTOOLS_X=1.13

LABSHR=/lab-share/IM-Gutierrez-e2/Public
SAMPLES=( $(cat ${LABSHR}/vitor/ase/bcell_bulk/bcell_samples.txt) )
BAMS=( "${SAMPLES[@]/#/$TEMP_WORK/gatk}" )
BAMS=( "${BAMS[@]/%/_recal.bam}" )
GENOME=${LABSHE}/vitor/ase/data/GRCh38.primary_assembly.genome.fa
OUTDIR=${TEMP_WORK}/gatk
PREFIX=${OUTDIR}/MG8989
MERGED=${PREFIX}.merged.bam

samtools merge -@ $SLURM_CPUS_PER_TASK $MERGED "${BAMS[@]}"

# HaplotypeCaller
VCF=${PREFIX}.vcf.gz

gatk HaplotypeCaller \
    -R $GENOME \
    -I $MERGED \
    -O $VCF \
    --dont-use-soft-clipped-bases true \
    -stand-call-conf 10 \	
    --native-pair-hmm-threads $SLURM_CPUS_PER_TASK

FILTVCF=${PREFIX}.filtered.vcf.gz
LOG=${PREFIX}.snpCluster.log 

gatk VariantFiltration \
    -R $GENOME \
    -V $VCF \
    -O $FILTVCF \
    -window 35 -cluster 3 &> $LOG 

# Save only SNVs
FILTVCF_SNV=${PREFIX}.filtered.biallelic.het.SNVs.vcf.gz 

bcftools norm -m -both $FILTVCF |\
    bcftools view --types snps - |\
    bcftools norm -m +snps - |\
    bcftools view -m2 -M2 --genotype ^hom -O z -o $FILTVCF_SNV

# Extract only SNVs reported by dbSNP
DBSNP=${LABSHR}/vitor/ase/data/dbsnp_155.hg38pri.vcf.gz
POS=${PREFIX}.pos
DBSNP_VCF=${PREFIX}.dbSNP.vcf.gz

# Extract dbSNP SNV variants for the current sample
zcat $FILTVCF_SNV |\
    grep -v "^#" |\
    awk '{ print $1"\t"$2 }' |\
    uniq > $POS

bcftools view -R $POS --types snps -O z -o $DBSNP_VCF $DBSNP

# Now refilter the sample VCF for the SNV variants present in dbSNP
FINAL_POS=${PREFIX}.final.pos

zcat $DBSNP_VCF |\
    grep -v "^#" |\
    awk '{ print $1"\t"$2 }' |\
    uniq > $FINAL_POS

FINAL_VCF=${LABSHR}/vitor/bcell_bulk/results/gatk/${SAMPLE}.knownSNVs.biallelic.het.vcf

tabix -p vcf $FILTVCF_SNV  
bcftools view -R $FINAL_POS -O v -o $FINAL_VCF $FILTVCF_SNV 

