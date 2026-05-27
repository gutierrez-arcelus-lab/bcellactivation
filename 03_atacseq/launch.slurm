#!/usr/bin/bash

#SBATCH --nodes=1                
#SBATCH --cpus-per-task=1        
#SBATCH --mem=4G                 
#SBATCH --time=72:00:00
#SBATCH --partition=bch-compute
#SBATCH --job-name=atac
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --output=atac_nfcore.%j.out
#SBATCH --error=atac_nfcore.%j.err

# load singularity
module purge
module load singularity

# set temp dirs
export NXF_SINGULARITY_CACHEDIR=${TEMP_WORK}/temp/singularity/nxf_cache
export SINGULARITY_TMPDIR=${TEMP_WORK}/temp/singularity/tmp
export NXF_TEMP=${TEMP_WORK}/temp/nextflow/tmp
export NXF_WORK=${TEMP_WORK}/temp/nextflow/work

# mv to execution dir
cd $SLURM_SUBMIT_DIR

# run
conda run -p $HOME/miniconda3/envs/nf-core nextflow run nf-core/atacseq -r 2.0 \
    -c e2.config \
    --input samplesheet.csv \
    --read_length 150 \
    --narrow_peak \
    --mito_name chrM \
    --blacklist /lab-share/IM-Gutierrez-e2/Public/References/encode/blacklist/hg38-blacklist.v3.bed \
    --outdir ./results \
    --multiqc_title B_cell_atac_nfcore \
    -resume

