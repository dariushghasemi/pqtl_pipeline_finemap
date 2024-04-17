suppressMessages(library(optparse))

option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
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

# separate path from seqid given by study_id, then make name pattern
sumstat_name <- basename(opt$study_id)
sumstat_path <- paste0(dirname(opt$study_id), "/")
file_pattern <- paste0(sumstat_name, "_chr(\\d+)_dataset_aligned.tsv.gz")

# Load-in summary statistics munged and aligned, binned for all chromosomes
#dataset_aligned_list <- list.files(pattern=paste0(opt$study_id, "_chr(\\d+)_dataset_aligned.tsv.gz"), path="./")
dataset_aligned_list <- list.files(pattern = file_pattern, path = sumstat_path) 

# add directory in which the files exist to the file names
dataset_aligned_list_full <- paste0(sumstat_path, dataset_aligned_list)

# Merge them
dataset_aligned <- as.data.frame(rbindlist(lapply(dataset_aligned_list_full, function(x){
  fread(x, data.table=F)
}))) %>% arrange(CHR)

# Save
#fwrite(dataset_aligned, paste0(opt$study_id, "_dataset_aligned.tsv.gz"), quote=F, na=NA, sep="\t")
fwrite(dataset_aligned, paste0(opt$outdir, "_dataset_aligned.tsv.gz"), quote=F, na=NA, sep="\t")


################
# Locus breaker
################

loci_list <- as.data.frame(rbindlist(
  lapply(dataset_aligned %>% group_split(phenotype_id), function(x){
    ### Check if there's any SNP at p-value lower than the set threshold. Otherwise stop here
    if(any(x %>% pull(p) < opt$p_thresh1)){
      ### Loci identification
      locus.breaker(
        x,
        p.sig=opt$p_thresh1,
        p.limit=opt$p_thresh2,
        hole.size=opt$hole,
        p.label="p",
        chr.label="CHR",
        pos.label="BP")
    }
  }) %>% discard(is.null)
))

### Add study ID to the loci table. Save
#loci_list <- loci_list %>% mutate(study_id=opt$study_id)
loci_list <- loci_list %>% mutate(study_id=sumstat_name)

#fwrite(loci_list, paste0(opt$study_id, "_loci.tsv"), sep="\t", quote=F, na=NA)
fwrite(loci_list, paste0(opt$outdir, "_loci.tsv"), sep="\t", quote=F, na=NA)

#cat(paste0("\n", nrow(loci_list), " significant loci identified for ", opt$study_id, "\n"))
cat(paste0("\n", nrow(loci_list), " significant loci identified for ", sumstat_name, "\n"))

cat("\nAnalysis is done!\n")


