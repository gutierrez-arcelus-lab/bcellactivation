#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=84,88,90
#SBATCH --job-name=vireo
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/vireo_%A_%a.out
#SBATCH -e /temp_work/ch229163/log/vireo_%A_%a.err

source /programs/biogrids.shrc
source ${HOME}/miniconda3/etc/profile.d/conda.sh

conda activate vireo

BATCH=19${SLURM_ARRAY_TASK_ID}
LABSHR=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq

if [ "$BATCH" -eq 1984 ]; then
    DIR=${LABSHR}/SN0268787/broad/hptmp/curtism/bwh10x/KW10456_Maria/221101_10X_KW10456_bcl/cellranger-6.1.1/GRCh38/BRI-${BATCH}/outs 
else
    DIR=${LABSHR}/SN0273471/broad/hptmp/sgurajal/bwh10x/KW10598_mgutierrez/221211_10X_KW10598_bcl/cellranger-6.1.1/GRCh38/BRI-${BATCH}_hashing/outs
fi

BAM=${DIR}/possorted_genome_bam.bam
BARCODE=${DIR}/filtered_feature_bc_matrix/barcodes.tsv.gz
VARIANTS=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/citeseq/demultiplexing/genotypes/data/allchroms_mgb_variants.vcf.gz
OUT=${SLURM_SUBMIT_DIR}/data/${BATCH}

mkdir -p $OUT

cellsnp-lite -s $BAM -b $BARCODE -R $VARIANTS -O $OUT -p 10 --minCOUNT 20 --gzip --genotype

VCF=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/citeseq/demultiplexing/genotypes/data/allchroms_mgb.vcf.gz
RES=${SLURM_SUBMIT_DIR}/results/${BATCH}

mkdir -p $RES

vireo -c $OUT -d $VCF -o $RES --genoTag=GT

conda deactivate
