#!/bin/bash

#SBATCH --job-name test_abf
#SBATCH --output %j_test_abf.log
#SBATCH --partition cpuq #-interactive
#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH --time 00:30:00

source ~/.bashrc
module -s load singularity/3.8.5

# run the pipeline
conda activate /exchange/healthds/software/envs/snakemake
snakemake   --profile slurm  #--unlock
conda deactivate
