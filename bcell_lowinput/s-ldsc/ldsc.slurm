#!/usr/bin/bash 

#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -J ldsc
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ldsc.%A-%a.out
#SBATCH -e /temp_work/ch229163/log/ldsc.%A-%a.err

source /programs/biogrids.shrc
source $HOME/miniconda3/etc/profile.d/conda.sh

conda activate ldsc

cd $SLURM_SUBMIT_DIR

GWAS=data/sumstats/PASS_Lupus.sumstats.gz
WEIGHTS=data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
FREQ=data/1000G_Phase3_frq/1000G.EUR.QC.

BASELINE=data/1000G_Phase3_baseline_v1.2_ldscores/baseline.
TURQ=data/gene_sets/ldscores/turquoise.
BLUE=data/gene_sets/ldscores/blue.
BROWN=data/gene_sets/ldscores/brown.
YELLOW=data/gene_sets/ldscores/yellow.
GREEN=data/gene_sets/ldscores/green.
RED=data/gene_sets/ldscores/red.
BLACK=data/gene_sets/ldscores/black.

python2.7 /lab-share/IM-Gutierrez-e2/Public/vitor/software/ldsc/ldsc.py \
    --h2 $GWAS \
    --ref-ld-chr $BASELINE,$TURQ,$BLUE,$BROWN,$YELLOW,$GREEN,$RED,$BLACK \
    --w-ld-chr $WEIGHTS \
    --frqfile-chr $FREQ \
    --print-coefficients \
    --out module_results

conda deactivate
