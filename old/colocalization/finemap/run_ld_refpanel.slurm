#!/usr/bin/sh

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=8:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-41
#SBATCH --job-name=LD1KGP
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/LD1KGP-%A-%a

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12
export R_X=4.1

IDS=${SLURM_SUBMIT_DIR}/data/europeans_samples.txt
REGIONSFILE=${SLURM_SUBMIT_DIR}/data/regions_bentham_1mb.tsv 
REGION=$( awk -v ARR="$SLURM_ARRAY_TASK_ID" 'FNR == ARR { print $1 }' $REGIONSFILE )
CHR=$( echo $REGION | cut -d':' -f1 )
VCF1K=/reference_databases/1000G_VCF/GRCh38/Phased_VCFs/CCDG_14151_B01_GRM_WGS_2020-08-05_${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz
PREFIX=${SLURM_SUBMIT_DIR}/data/${REGION}
VCFOUT=${PREFIX}.vcf.gz
VCFADJ=${PREFIX}.adj.vcf.gz
VARIANTS=${PREFIX}.variants.txt

bcftools view --threads $SLURM_CPUS_PER_TASK -r $REGION $VCF1K |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --samples-file $IDS --force-samples |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --min-ac 1:minor - |\
    bcftools view --threads $SLURM_CPUS_PER_TASK -e INFO/AC==INFO/AN - |\
    bcftools norm --threads $SLURM_CPUS_PER_TASK -m +both - |\
    bcftools annotate --threads $SLURM_CPUS_PER_TASK -x FORMAT - |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --genotype ^miss -O z -o $VCFOUT -

Rscript ${SLURM_SUBMIT_DIR}/adjust_vcf.R $VCFOUT $VARIANTS $VCFADJ

plink --vcf $VCFADJ --r square spaces --keep-allele-order --out ${PREFIX}
rm ${PREFIX}.nosex
