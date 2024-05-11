

gwas=/exchange/healthds/pQTL/results/INTERVAL/chunk_7023/chunk_7023_output/chunk_7023/results/gwas/seq.9832.33.gwas.regenie.gz

#tabix -h $gwas 22:43824730-44824730 | locuszoom --metal - markercol ID --epacts-beg-col GENPOS --epacts-end-col GENPOS --epacts-chr-col CHROM --refsnp chr22:44324730  --chr 22  --start  43824730  --end 44824730  --build hg19 --no-ld -plotonly --prefix "11-Mar-24_22_44324730"

locuszoom  \
    --metal $gwas  \
    --markercol ID \
    --epacts-beg-col GENPOS \
    --epacts-end-col GENPOS \
    --epacts-chr-col CHROM \
    --refsnp chr22:44324730  \
    --chr 22  \
    --flank 5kbp  \
    --pop EUR \
    --build hg19 \
    --source 1000G_March2012 \
    --plotonly   \
    --prefix "11-Mar-24_22_44324730"