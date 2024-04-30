FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="e5e15a739a807ed5368d20504462c433a3036490881814232e4aee2bf799dda1"

# Step 1: Retrieve conda environments


# Conda environment:
#   source: workflow/envs/r_environment.yml
#   prefix: /conda-envs/baff2ac264f1f12f7171aed64ac3702e
#   name: r_finemap
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#     - R
#   dependencies: #R>=4.3
#     - r-base=4.3.3
#     - r-optparse=1.7.4
#     - r-tidyverse=2.0.0
#     - r-R.utils=2.12.3
#     - r-coloc=5.2.3
#     - r-cowplot=1.1.3
#     - r-corrplot=0.92
#     - r-bigsnpr=1.12.2
#     - r-stringi=1.8.3
#     - r-patchwork=1.2.0
#     - r-plyr=1.8.9
#     - r-reshape2=1.4.4
#     - r-RColorBrewer=1.1-3
#     - r-igraph=2.0.2
#     - r-Matrix=1.6-5
#     - bioconductor-Gviz=1.46.1
#     - bioconductor-EnsDb.Hsapiens.v75=2.99.0
RUN mkdir -p /conda-envs/baff2ac264f1f12f7171aed64ac3702e
COPY workflow/envs/r_environment.yml /conda-envs/baff2ac264f1f12f7171aed64ac3702e/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/baff2ac264f1f12f7171aed64ac3702e --file /conda-envs/baff2ac264f1f12f7171aed64ac3702e/environment.yaml && \
    mamba clean --all -y
