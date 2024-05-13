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

## Regional association plot of leading pQTLs

- To use LocusZoom standalone tool v1.4, we needed to convert the genotype file from pgen/pvar/psam to VCF format. We used PLINK v2.00a5LM to do the conversion (Sat, 21:57, 11-May-24).

``` bash
# convert pgen to bg-zipped vcf
plink2 --pfile /scratch/giulia.pontali/test_snakemake/genomics_QC_pipeline/results/pgen/impute_recoded_selected_sample_filter_hq_var_19 --recode vcf bgz --out interval.imputed.info70.chr22

# create an index file
tabix -p vcf 

# modify the summary results for locuszoom
zcat seq.9832.33.gwas.regenie.gz | head -1 > my_header.txt
tabix seq.9832.33.gwas.regenie.gz 22:44324530-44324930 | sed -E 's/([0-9]+:[0-9]+):([A-Z]):([A-Z])/\1_\2\/\3/g' > my_region.txt

# all together
cat my_header.txt > chr22.txt && cat my_region.txt >> chr22.txt && cat chr22.txt | bgzip > chr22.txt.gz && tabix chr22.txt.gz -s1 -b2 -e2 -f

# subset vcf
bcftools view interval.imputed.info70.chr22.vcf.gz  -r 22:44324530-44324930 -Oz -o interval.imputed.info70.chr22.cut.vcf.gz

# subset vcf via tabix
tabix -p vcf interval.imputed.info70.chr22.vcf.gz 22:43824730-44824730 --print-header | sed -E 's/([0-9]+:[0-9]+):([A-Z]):([A-Z])/\1_\2\/\3/g' | bgzip > interval.imputed.info70.chr22.cut.aligned.vcf.gz

tabix -p vcf interval.imputed.info70.chr22.cut.aligned.vcf.gz

# change colon with underscore in vcf file 22:44324530-44324930
#bcftools view interval.imputed.info70.chr22.cut.vcf.gz | sed -E 's/([0-9]+:[0-9]+):([A-Z]):([A-Z])/\1_\2\/\3/g' | bgzip > interval.imputed.info70.chr22.cut.aligned.vcf.gz

# save snps list in vcf file
tabix interval.imputed.info70.chr22.cut.aligned.vcf.gz 22:43824730-44824730 | cut -f3 > snp.list

# compute ld
plink --vcf interval.imputed.info70.chr22.cut.aligned.vcf.gz       --ld-snp-list snp.list        --ld-window 10000    --ld-window-kb 250    --r2 dprime        --ld-window-r2 0        --out ld_region

# compute ld with lead variant
plink --vcf interval.imputed.info70.chr22.cut.aligned.vcf.gz       --ld-snp 22:44324730_T/C      --ld-window 10000        --ld-window-kb 500        --r2 dprime        --ld-window-r2 0        --out ld_region

# reform ld file for LZ
cat ld_region.ld | awk '{OFS="\t"; {print $3, $6, $8, $7}}' | sed -e '1s/SNP_A/snp1/' -e '1s/SNP_B/snp2/' -e '1s/DP/dprime/' -e '1s/R2/rsquare/' > my_ld.txt

# locuszoom with LD from VCF
tabix chr22.txt.gz 22:44324530-44324930 -h --print-header | locuszoom --metal - --markercol ID --pvalcol LOG10P --no-transform --refsnp 22:44324730  --flank 200  --build hg19 --ld-vcf interval.imputed.info70.chr22.cut.aligned.vcf.gz  --plotonly --prefix "13-Mar-24_ld_vcf"

# locuszoom with precalculated LD and via STDIN *** working ***
tabix chr22.txt.gz 22:43824530-44824930 -h --print-header | locuszoom --metal - --markercol ID --pvalcol LOG10P --no-transform --refsnp 22:44324730  --flank 250kbp  --build hg19 --ld my_ld.txt  --ld-measure rsquare --plotonly --prefix "13-Mar-24_ld_user"

```

- A regional association plot for an example locus '22:44324730' associated with 'seq.9832.33' protein was generated using a user-defined local LD in Interval study and save [here](/home/dariush.ghasemi/projects/pqtl_pipeline_finemap/lz_plot/13-Mar-24_ld_user_240513_22_44324730.pdf) (Mon, 19:30, 13-May-24).

- Successfully tested running of cojo after merging with mapping via pipeline on clusters (Tue, 12:50, 14-May-24).

Dariush