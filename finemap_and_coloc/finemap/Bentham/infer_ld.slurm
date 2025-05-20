#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=8:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH --array=1-36
#SBATCH --job-name=kgp-ld
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/kgpld-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/kgpld-%A-%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.20

CPUS=$SLURM_CPUS_PER_TASK

cd ${SLURM_SUBMIT_DIR}

# Window
REGION=$( awk -v ARR="$SLURM_ARRAY_TASK_ID" 'FNR == ARR { print $2 }' ./data/windows.tsv )
CHR=$( echo $REGION | cut -d':' -f1 )

IDS=./data/ref_panel.txt

if [[ "$CHR" == X ]]; then
    VCF=/reference_databases/1000G_VCF/phase3/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz
else
    VCF=/reference_databases/1000G_VCF/phase3/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
fi

# OUT
PREFIX=./data/ld/window${SLURM_ARRAY_TASK_ID}
VCFOUT=${PREFIX}.vcf.gz

# VCF processing
bcftools view --threads $CPUS -r $REGION $VCF |\
    bcftools view --threads $CPUS --samples-file $IDS --force-samples - |\
    bcftools norm --threads $CPUS -m - - |\
    bcftools view --threads $CPUS --min-ac 2:minor - |\
    bcftools view --threads $CPUS --min-ac 2:major - |\
    bcftools annotate --threads $CPUS -x QUAL,FILTER,INFO,FORMAT -O z -o $VCFOUT -

# LD
plink --vcf $VCFOUT --r square spaces --real-ref-alleles --threads $CPUS --out ${PREFIX}_r
rm ${PREFIX}_r.nosex
