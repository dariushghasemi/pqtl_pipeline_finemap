# pqtl_pipeline_finemap
Fine mapping analysis within the pQTL pipeline project

```{bash}

# get to the project dir
cd projects/pqtl_pipeline_finemap/

# copy R scripts into home directory
cp /exchange/healthds/pQTL/pQTL_workplace/fine_mapping_genomics/nf-hcoloc/bin/*.R projects/pqtl_pipeline_finemap/workflow/scripts/

# create conda env to install and load R libraries
mkdir workflow/envs
cd envs/
touch r_finemap.yml

# add packages with their version in yaml file
conda env create -f r_environment.yml 

# test r code
conda activate r_finemap
Rscript ../scripts/s01_sumstat_munging.R 
conda deactivate

# get to the env to install missed libraries
cd ../envs/
conda env update -f r_environment.yml 
conda env list

# try to reactivate env to test it
conda activate r_finemap
Rscript ../scripts/s01_sumstat_munging.R
conda deactivate

```

- First step is tested. Second is ongoing! (19:00, Fri, 12-Apr-24)

- Configuration and snakefile were created to integrate the codes into the pipeline (23:50, Sun, 14-Apr-24)!

- Lunching the three first steps of the workflow on slurm integrated into snakemake pipeline (03:15, Mon, 15-Apr-24).

- Step four was also debugged after integrating t

Dariush