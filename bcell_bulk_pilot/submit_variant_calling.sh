#!/bin/bash

JOB1=$( sbatch --parsable star_1st.slurm )

JOB2=$( sbatch --parsable --dependency=afterok:$JOB1 star_2nd.slurm )

JOB3=$( sbatch --parsable --dependency=afterok:$JOB2 gatk.slurm )

JOB4=$( sbatch --parsable --dependency=afterok:$JOB3 haplotypecaller.slurm )

./concat_vcf.sh

JOB5=$( sbatch --parsable --dependency=afterok:$JOB4 star_wasp.slurm )

JOB6=$( sbatch --parsable --dependency=afterok:$JOB5 asereadcounter.slurm )
