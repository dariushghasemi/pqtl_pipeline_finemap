suppressMessages(library(optparse))

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--input", default=NULL, help="Path and file name of GWAS summary statistics"),
  make_option("--is_molQTL", default=NULL, help="Whether the summary statistics provided are for eQTLs, ATAC, ChipSeq etc. data"),
  make_option("--key", default=NULL, help="If GWAS is from a molQTL, name of column reporting the trait/phenoptype"),
  make_option("--chr_lab", default="CHROM", help="Name of chromosome column in GWAS summary statistics"),
  make_option("--pos_lab", default="GENPOS", help="Name of genetic position column in GWAS summary statistics"),
  make_option("--rsid_lab", default="ID", help="Name of rsid column"),  
  make_option("--a1_lab", default="ALLELE1", help="Name of effect allele column"),  
  make_option("--a0_lab", default="ALLELE0", help="Name of NON effect allele column"),  
  make_option("--freq_lab", default="A1FREQ", help="Name of effect allele frequency column"),  
  make_option("--n_lab", default="N", help="Name of sample size column"),  
  make_option("--effect_lab", default="BETA", help="Name of effect size column"),  
  make_option("--se_lab", default="SE", help="Name of standard error of effect column"), 
  make_option("--pvalue_lab", default="P", help="Name of p-value of effect column"),
  make_option("--type", default=NULL, help="Type of phenotype analysed - either 'quant' or 'cc' to denote quantitative or case-control"),
  make_option("--sdY", default=NULL, help="For a quantitative trait (type==quant), the population standard deviation of the trait. If not given, it will be estimated beta and MAF"),
  make_option("--s", default=NULL, help="For a case control study (type==cc), the proportion of samples in dataset 1 that are cases"),
  make_option("--study_id", default=NULL, help="Id of the study")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))

## Throw error message - GWAS summary statistics file MUST be provided!
if(is.null(opt$input)){
  print_help(opt_parser)
  stop("Please specify the path and file name of your GWAS summary statistics in --path option", call.=FALSE)
}

##################################
# Load in and munge GWAS sum stats
##################################
gwas <- fread(opt$input, data.table = F)

# If the trait column is NOT provided, add the same one for the whole sum stat
if(opt$is_molQTL==FALSE){
  gwas <- gwas %>% mutate(phenotype_id="full")
  opt$key="phenotype_id"
}

# If the trait column is provided, simply rename it
if(opt$is_molQTL==TRUE){
  gwas <- gwas %>% rename(phenotype_id=opt$key)
}

dataset <- dataset.munge_hor(
  gwas
  ,snp.lab = opt$rsid_lab
  ,chr.lab = opt$chr_lab
  ,pos.lab = opt$pos_lab
  ,a1.lab = opt$a1_lab
  ,a0.lab = opt$a0_lab
  ,beta.lab = opt$effect_lab
  ,se.lab = opt$se_lab
  ,freq.lab = opt$freq_lab
  ,pval.lab = opt$pvalue_lab
  ,n.lab = opt$n_lab
  ,type = opt$type
  ,sdY = opt$sdY
  ,s = opt$s
)

#### Remove multi-allelic SNPs in the same position! All occurrences
dataset <- as.data.frame(dataset %>% group_by(phenotype_id, CHR, BP) %>% filter(n() == 1))

saveRDS(dataset, file=paste0(opt$study_id, "_dataset.rds"))
cat(paste0("\nGWAS summary statistics for ", opt$study_id, " has been harmonised!\n"))



