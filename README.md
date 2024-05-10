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

- The first two steps of the original workflow are not needed for the pQTLs pipeline. Therefore, it is better to start with the step 3, locus breaker, directly using protein GWAS summary stats from REGENIE which is already aligned (Fri, 11:58, 26-Apr-24). 

- Locus breaker is incorporated into the pipeline. Time to do so with COJO execution (Sat, 03:45, 27-Apr-24).

- After that the data manager prepared the pipeline for execution on the proteins sequence list, I created a new branch `git checkout -b develop` and switched to it for further adjusting the pipeline for COJO execution (Tue, 13:18, 30-Apr-24).

- The rule of COJO-ABF fine-mapping was modified properly to work with GWAS results including -LOG10P. The rule gets executed after a successful run of the locus breaker rule in the pipeline. The changes in the development branch need to be merged with the main branch later. Before lunching pipeline, we need to incorporate mapping file in run_COJO rule to tackle the issue related to GWAS results file with alleles in alphabetical order (Wed, 11:50, 01-May-24).

- The list of proteins in order to test the pipeline (Fri, 14:45, 03-May-24).

```bash
# take path of instance proteins and list them in config folder of the pipeline 
find /exchange/healthds/pQTL/results/INTERVAL/*/*/*/results/gwas/*.gz | grep -E '19819.7|12708.91|12730.3|7930.3|7935.26|7943.16|13124.20' > path_list.txt /scratch/.../projects/conf/path_list.txt

# full list of the proteins for final run
ls /exchange/healthds/pQTL/results/INTERVAL/chunk_*/chunk_*_output/chunk_*/results/gwas/seq.*.gwas.regenie.gz > /scratch/.../projects/conf/path_list.txt
```

- Locus breker results were integrated and sent to the PI (20:06, Mon, 06-May-24).
- After screening the results of locus breaker function, it turned out that there are many loci whose width were so lengthy or even equal to zero. Need to check its performance once applied to the meta-analysis and then draw LocusZoom plots (Fri, 17:30, 10-May-24).  

Dariush