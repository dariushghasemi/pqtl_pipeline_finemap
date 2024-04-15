#!/usr/bin/bash

#SBATCH --job-name worker
#SBATCH --output %j_test_s03.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1
#SBATCH --mem 8G
#SBATCH --time 8:00:00

#    --input test_pwas.regenie.tsv.gz  \
pwas=/exchange/healthds/pQTL/results/INTERVAL_OLD/harmonized_sumstats/outputs/seq.7930.3/seq.7930.3.regenie.tsv.gz

#RAM needed: 16GB

# Rscript s01_sumstat_munging.R  \
#     --pipeline_path /home/dariush.ghasemi/projects/pqtl_pipeline_finemap/workflow/scripts/  \
#     --input /exchange/healthds/pQTL/results/INTERVAL_OLD/harmonized_sumstats/outputs/seq.7930.3/seq.7930.3.regenie.tsv.gz  \
#     --is_molQTL FALSE  \
#     --type quant  \
#     --s NA  \
#     --study_id seq.7930.3


# Rscript s02_sumstat_alignment.R  \
#     --pipeline_path /home/dariush.ghasemi/projects/pqtl_pipeline_finemap/workflow/scripts/  \
#     --dataset seq.7930.3_dataset.rds  \
#     --study_id seq.7930.3  \
#     --grch 38  \
#     --chr_tabix 22


# Rscript s03_locus_breaker.R  \
#     --pipeline_path /home/dariush.ghasemi/projects/pqtl_pipeline_finemap/workflow/scripts/  \
#     --study_id seq.7930.3  \
#     --p_thresh1 5e-06  \
#     --p_thresh2 1e-04

Rscript scripts/s03_locus_breaker.R  \
    --pipeline_path /home/dariush.ghasemi/projects/pqtl_pipeline_finemap/workflow/scripts/ \
    --p_thresh1 5e-06 --p_thresh2 1e-04  \
    --study_id /home/dariush.ghasemi/projects/pqtl_pipeline_finemap/output/align/seq.13124.20  \
    --outdir /home/dariush.ghasemi/projects/pqtl_pipeline_finemap/output/break/seq.13124.20/seq.13124.20