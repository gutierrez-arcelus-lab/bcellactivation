#!/usr/bin/sh

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-39
#SBATCH --job-name=LD1KGP
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/LD1KGP-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/LD1KGP-%A-%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.20
export R_X=4.1

cd $SLURM_SUBMIT_DIR

# I/O
IDS=./data/eur_kgp.txt
REGIONLIST=data/langefeld_regions.tsv
REGION=$( awk -v ARR="$SLURM_ARRAY_TASK_ID" 'FNR == ARR { print $1 }' $REGIONLIST )
CHR=$( echo $REGION | cut -d':' -f1 )
VCF1K=/reference_databases/1000G_VCF/phase3/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
PREFIX=data/chr${REGION}
VCFOUT=${PREFIX}.vcf.gz

bcftools view --threads 4 -r $REGION $VCF1K |\
    bcftools view --threads 4 --samples-file $IDS --force-samples |\
    bcftools norm --threads 4 -m - - |\
    bcftools view --threads 4 --min-ac 2:minor - |\
    bcftools view --threads 4 --min-ac 2:major - |\
    bcftools annotate --threads 4 -x QUAL,FILTER,INFO,FORMAT -O z -o $VCFOUT -

plink --vcf $VCFOUT --r square spaces --real-ref-alleles --threads 4 --out ${PREFIX}_r
plink --vcf $VCFOUT --r2 square spaces --real-ref-alleles --threads 4 --out ${PREFIX}_r2
