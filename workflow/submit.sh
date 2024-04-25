#!/bin/bash

#SBATCH --job-name loop_abf
#SBATCH --output %j_loop_abf.log
#SBATCH --partition cpuq #-interactive
#SBATCH --cpus-per-task 1
#SBATCH --mem 4G
#SBATCH --time 10-00:00:00

source ~/.bashrc
module -s load singularity/3.8.5

# run the pipeline
conda activate /exchange/healthds/software/envs/snakemake
snakemake   --profile slurm  #--unlock
conda deactivate
