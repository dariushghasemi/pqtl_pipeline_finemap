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

- Step four was also debugged after integrating to the pipeline (16:30, Mon, 22-Apr-24).

- It was critical and needed to loop the fine-mapping analysis (step 4) over the identified index variants (P<5e-8, step 3 or locus breaker output) for each protein sequence. This action was done inside rule 'run_cojo' which allowed iterating the R script `s04_cojo_finemapping.R` inside the rule rather than looping the rule through the index variants. Therefore, we read TSV file containing the index variants and iterating the R script only for those variants. So that if there is no genome-wide significant variants for a protein sequence, rule 'run_cojo' would fail for that proetin.

- According to the abovementioned explanations, fine-mapping outputs including independet variants, credible set, regional association plots were saved for 3 out of 6 tested protein sequence in `~/projects/pqtl_pipeline_finemap/output/cojo/<seqid>` (Thu, 23:45, 25-Apr-24).

Dariush