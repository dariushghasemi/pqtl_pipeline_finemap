suppressMessages(library(optparse))

option_list <- list(
  make_option("--pipeline_path", default=NULL, help="Path where Rscript lives"),
  make_option("--dataset", default=NULL, help="Munged dataset in rds format"),
  make_option("--chr_tabix", default=NULL, help="Harmonization by tabix will be performed by chromosome "),
  make_option("--grch", default=NULL, help="Genomic build of GWAS summary statistics"),
  make_option("--tabix_bin", default="/ssu/gassu/software/htslib-tools/1.14/tabix", help="Path to tabix software"),
  make_option("--study_id", default=NULL, help="Id of the study")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Source function R functions
source(paste0(opt$pipeline_path, "funs_locus_breaker_cojo_finemap_all_at_once.R"))

cat("Genomic build is:", opt$grch, "\n\n")

#########################
# Harmonisation with map!
#########################

# Set universal mapping files (created by Edo, 1000genomes+HRC+TOPMed)
#if(opt$grch==38){ mappa="https://zenodo.org/records/10142577/files/Full_variant_map.GRCh38_sorted.tsv.gz" }
#if(opt$grch==37){ mappa="https://zenodo.org/records/10142577/files/Full_variant_map.GRCh37_sorted.tsv.gz" }

if(opt$grch==38){ mappa="/ssu/bsssu/ghrc38_reference/Full_variant_map.GRCh38_sorted.tsv.gz" }
if(opt$grch==37){ mappa="/ssu/bsssu/ghrc37_reference/Full_variant_map.GRCh37_sorted.tsv.gz" }

dataset_munged <- readRDS(opt$dataset)

# First check if opt$chr_tabix value is present in the summary statistics - otherwise skip the whole process
if(opt$chr_tabix %in% unique(dataset_munged$CHR)){

  dataset_aligned_tmp <- dataset.align(
    dataset_munged,
    study_id=opt$study_id,
    mappa,
    chr_tabix=opt$chr_tabix,
    tabix_bin=opt$tabix_bin,
    grch=opt$grch
  )
  gc()

  fwrite(dataset_aligned_tmp, paste0(opt$study_id, "_chr", opt$chr_tabix,"_dataset_aligned.tsv.gz"), sep="\t", quote=F, na=NA)

  cat(paste0("\nGWAS summary statistics for ", opt$study_id, "in chromosome ", opt$chr_tabix," have been aligned\n"))
}
