FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="554edbb710da5a2b8cbb45a61e2f586edded22fb58489eed03ac72773f61099a"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: envs/r_environment.yml
#   prefix: /conda-envs/50a3fe91d9589eca4a256154c89f39b6
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
RUN mkdir -p /conda-envs/50a3fe91d9589eca4a256154c89f39b6
COPY envs/r_environment.yml /conda-envs/50a3fe91d9589eca4a256154c89f39b6/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/50a3fe91d9589eca4a256154c89f39b6 --file /conda-envs/50a3fe91d9589eca4a256154c89f39b6/environment.yaml && \
    mamba clean --all -y
