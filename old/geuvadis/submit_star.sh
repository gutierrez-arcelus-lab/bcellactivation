#! /bin/bash

JOB1=$( sbatch --parsable ./star_1st.slurm )

sbatch --dependency=afterok:$JOB1 ./star_2nd.slurm
