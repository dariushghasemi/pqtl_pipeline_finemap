suppressMessages(library(optparse))

option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--input", default=NULL, help="Path and file name of GWAS summary statistics"),
  make_option("--p_thresh1", default=5e-08, help="Significant p-value threshold for top hits"),
  make_option("--p_thresh2", default=1e-05, help="P-value threshold for loci borders"),  
  make_option("--hole", default=250000, help="Minimum pair-base distance between SNPs in different loci"),
  make_option("--study_id", default=NULL, help="Id of the study"),
  make_option("--outdir", default=NULL, help="Output directory")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))
cat("\n\nChecking input file...")
## Throw error message - GWAS summary statistics file MUST be provided!
if(is.null(opt$input)){
  print_help(opt_parser)
  stop("Please specify the path and file name of your GWAS summary statistics in --path option", call.=FALSE)
}

##################################
# Load in and munge GWAS sum stats
##################################

cat("\n\nReading GWAS file...")
gwas <- data.table::fread(opt$input, data.table = F)
cat("done!\n")
# separate path from seqid given by study_id, then make name pattern
# sumstat_name <- basename(opt$study_id)
# sumstat_path <- paste0(dirname(opt$study_id), "/")
# file_pattern <- paste0(sumstat_name, "_chr(\\d+)_dataset_aligned.tsv.gz")

# Load-in summary statistics munged and aligned, binned for all chromosomes
#dataset_aligned_list <- list.files(pattern=paste0(opt$study_id, "_chr(\\d+)_dataset_aligned.tsv.gz"), path="./")
##dataset_aligned_list <- list.files(pattern = file_pattern, path = sumstat_path) 

# add directory in which the files exist to the file names
##dataset_aligned_list_full <- paste0(sumstat_path, dataset_aligned_list)

# Merge them
# dataset_aligned <- as.data.frame(rbindlist(lapply(dataset_aligned_list_full, function(x){
#   fread(x, data.table=F)
# }))) %>% arrange(CHR)

# Save
#fwrite(dataset_aligned, paste0(opt$study_id, "_dataset_aligned.tsv.gz"), quote=F, na=NA, sep="\t")
#fwrite(dataset_aligned, paste0(opt$outdir, "_dataset_aligned.tsv.gz"), quote=F, na=NA, sep="\t")


################
# Locus breaker
################

cat("\n\nDefining locus...")

# loci_list <- as.data.frame(rbindlist(
#   lapply(dataset_aligned %>% group_split(phenotype_id), function(x){
#     ### Check if there's any SNP at p-value lower than the set threshold. Otherwise stop here
#     if(any(x %>% pull(p) < opt$p_thresh1)){
#       ### Loci identification
#       locus.breaker(
#         x,
#         p.sig=opt$p_thresh1,
#         p.limit=opt$p_thresh2,
#         hole.size=opt$hole,
#         p.label="p",
#         chr.label="CHR",
#         pos.label="BP")
#     }
#   }) %>% discard(is.null)
# ))

# function to find the index variants at each locus
check_signif <- function(x){
  ### Check if there's any SNP at p-value lower than the set threshold. Otherwise stop here
  if(any(x %>% pull(LOG10P) > -log10(opt$p_thresh1))){
  ### Loci identification
  locus.breaker(
    x,
    p.sig     = -log10(opt$p_thresh1),
    p.limit   = -log10(opt$p_thresh2),
    hole.size = opt$hole,
    p.label   = "LOG10P",
    chr.label = "CHROM", 
    pos.label = "GENPOS"
  )
  }
} %>% discard(is.null)

#list of index variants
loci_list <- check_signif(gwas)
cat("done!\n")

### Add study ID to the loci table. Save
#loci_list <- loci_list %>% mutate(study_id=opt$study_id)
loci_list$study_id <- basename(opt$study_id)

cat("\n\nSaving index variants...")
#fwrite(loci_list, paste0(opt$study_id, "_loci.tsv"), sep="\t", quote=F, na=NA)
fwrite(loci_list, paste0(opt$outdir, "_loci.tsv"), sep="\t", quote=F, na=NA)
cat("done!\n")

#cat(paste0("\n", nrow(loci_list), " significant loci identified for ", opt$study_id, "\n"))
cat(paste0("\n", nrow(loci_list), " significant loci identified for ", sumstat_name, "\n"))

cat("\n\nLocus breaker is done!\n\n")


