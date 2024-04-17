suppressMessages(library(optparse))
#scipen=0

# Get arguments specified in the sbatch
option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--chr", default=NULL, help="Locus chromosome"),
  make_option("--start", default=NULL, help="Locus starting position"),
  make_option("--end", default=NULL, help="Locus ending position"),
  make_option("--phenotype_id", default=NULL, help="Trait for which the locus boundaries have been identified - relevant in cases of molQTLs"),
  make_option("--dataset_aligned", default=NULL, help="GENOME-WIDE munged and aligned dataset file"),
  make_option("--p_thresh3", default=1e-04, help="Noise p-values threshold for COJO"),
  make_option("--maf", default=1e-04, help="MAF filter", metavar="character"),
  make_option("--bfile", default=NULL, help="Path and prefix name of custom LD bfiles (PLINK format .bed .bim .fam)"),
  make_option("--plink2_bin", default="/ssu/gassu/software/plink/2.00_20211217/plink2", help="Path to plink2 software"),
  make_option("--gcta_bin", default="/ssu/gassu/software/GCTA/1.94.0beta/gcta64", help="Path to GCTA software"),
  make_option("--p_thresh4", default=1e-06, help="P-value significant threshold for redefining loci boundaries post-COJO"),  
  make_option("--hole", default=250000, help="Minimum pair-base distance between SNPs in different loci"),
  make_option("--cs_thresh", default=99, help="Percentage of credible set"),
  make_option("--study_id", default=NULL, help="Id of the study"),
  make_option("--nf_hcoloc_v", default=NULL, help="Version of nf hcoloc pipleine used, for reporting sake")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))

# Slightly enlarge locus by 200kb!
locus_name <- paste0(opt$chr, "_", opt$start, "_", opt$end)
opt$chr <- as.numeric(opt$chr)
opt$start <- as.numeric(opt$start) -100000
opt$end <- as.numeric(opt$end) + 100000

# GWAS input
dataset_aligned <- fread(opt$dataset_aligned, data.table=F) %>% dplyr::filter(phenotype_id==opt$phenotype_id)


###############
# Perform cojo
###############

# Break dataframe in list of single rows
conditional.dataset <- cojo.ht(
  D=dataset_aligned,
  locus_chr=opt$chr,
  locus_start=opt$start,
  locus_end=opt$end,
  p.thresh=as.numeric(opt$p_thresh3),
  maf.thresh=as.numeric(opt$maf),
  bfile=opt$bfile,
  gcta.bin=opt$gcta_bin,
  plink.bin=opt$plink2_bin
)

# Plot conditioned GWAS sum stats
pdf(paste0(opt$study_id, "_locus_chr", locus_name, "_conditioned_loci.pdf"), height=3.5*nrow(conditional.dataset$ind.snps), width=10) ### have the original loci boundaries in the name, or the slightly enlarged ones?
plot.cojo.ht(conditional.dataset) + plot_annotation(paste("Locus chr", locus_name))
dev.off()



####################
# Locus breaker BIS
###################

### Repeat only on dataset that have been conditioned!!
conditional.dataset$results <- lapply(conditional.dataset$results, function(x){
  
  ### Check if there's any SNP at p-value lower than the set threshold. Otherwise stop here
  if(isTRUE(any(x %>% pull(pC) < opt$p_thresh4))){
    new_bounds <- locus.breaker(
      x,
      p.sig=as.numeric(opt$p_thresh4),
      p.limit=as.numeric(opt$p_thresh3),
      hole.size=opt$hole,
      p.label="pC",
      chr.label="Chr",
      pos.label="bp")
    
    # Slightly enlarge locus by 200kb!
    new_bounds <- new_bounds %>% dplyr::mutate(start=as.numeric(start)-100000, end=as.numeric(end)+100000)
    
    # Remove SNPs not included in loci boundaries
    x %>% filter(bp >= new_bounds$start & bp <= new_bounds$end)
  }
})

## Remove eventually empty dataframes (caused by p_thresh4 filter)  
conditional.dataset$results <- conditional.dataset$results %>% discard(is.null)




#############
# Finemapping
#############

# Perform finemapping of each conditional dataset
finemap.res <- lapply(conditional.dataset$results, function(x){
  finemap.cojo(x, cs_threshold=opt$cs_thresh)
})


#########################################
# Organise list of what needs to be saved
#########################################

lapply(finemap.res, function(x){
  file_name <- paste0(opt$study_id, "_", opt$phenotype_id, "_locus_chr", locus_name, "_top_snp_", unique(x$cojo_snp), "_finemap")
  if(opt$phenotype_id=="full"){file_name <- gsub("_full", "", file_name)}
  
  # .rds object collecting 1) lABF, 2) beta, 3) pos for all SNPs, 3) list of SNPs in the credible set
  saveRDS(x %>% select(-cojo_snp), file=paste0(file_name, ".rds"))
  
  # .tsv with 1) study id and trait (if molQTL) locus info, 2) list of SNPs in the 99% credible set and 3) path and name of correspondent .rds file --> append each row to a master table collecting all info from processed sum stats
  ### Idea: create guidelines for study ids
  tmp <- data.frame(
    study_id=opt$study_id,
    phenotype_id=ifelse(opt$phenotype_id=="full", NA, opt$phenotype_id),
    credible_set=paste0(x %>% filter(is_cs==TRUE) %>% pull(snp), collapse=","),
    path=paste0(
      gsub("(.*)/work/.*", "\\1", getwd()), #### Nextflow working directory "work" hard coded - KEEP in mind!!
      "/results/finemap/", file_name, ".rds"),
    nf_hcoloc_v=opt$nf_hcoloc_v  
  )
  fwrite(tmp, paste0(file_name, "_coloc_info_table.tsv"), sep="\t", quote=F, col.names = F, na=NA)
})
