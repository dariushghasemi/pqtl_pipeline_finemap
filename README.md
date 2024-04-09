# pqtl_pipeline_finemap
Fine mapping analysis within the pQTL pipeline project

```bash
cp /exchange/healthds/pQTL/pQTL_workplace/fine_mapping_genomics/nf-hcoloc/bin/*.R projects/pqtl_pipeline_finemap/
  962  cd projects/pqtl_pipeline_finemap/
  968  mkdir workflow
  969  mkdir workflow/scripts
  970  mv s0*.R workflow/scripts/
  971  cd workflow/scripts/
  976  Rscript scripts/s01_sumstat_munging.R 
  977  R scripts/s01_sumstat_munging.R 
  978  ./scripts/s01_sumstat_munging.R 
  979  mkdir workflow/envs
  980  mkdir envs
  983  cd envs/
  984  touch r_finemap.yml
  995  conda env create -f r_environment.yml 
  996  exit
  997  cd projects/pqtl_pipeline_finemap/workflow/envs/
  998  conda env create -f r_environment.yml 
  999  conda activate r_finemap
 1002  cd ../scripts/
 1004  ./s01_sumstat_munging.R 
 1005  Rscript s01_sumstat_munging.R 
 1006  conda deactivate
 1007  cd ../envs/
 1008  conda env update -f r_environment.yml 
 1009  conda env list
 1010  conda activate r_finemap
 1011  Rscript s01_sumstat_munging.R 
 1012  cd ../scripts/
 1013  Rscript s01_sumstat_munging.R 
 1014  conda deactivate


```