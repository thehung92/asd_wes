#!/usr/bin/env Rscript
#
library(data.table)
library(tidyverse)
# read sample name from vcf file
INPUT="temp/data/asd.hg38_rename.vcf.gz"
QUERY="bcftools query -l"
CMD=paste(QUERY, INPUT)
vt.sample <- system(CMD, intern = TRUE)
# read pedigree data from IBD-clean file
INPUT="./output/asd.hg38_ibd-clean.fam"
fam <- fread(file=INPUT, header=FALSE,
             col.names = c("FID", "IID", "PID", "MID", "SEX", "PHENO"))
# exclude family of 2 members because these will not be helpful in denovo detection
vt.exclude <- count(fam, FID) %>% filter(n==2) %>% .$FID
fam <- fam %>%
  filter(! FID %in% vt.exclude)
# rename sample based on vcf sample name
fam <- fam %>%
  mutate_at(c("IID","PID","MID"), ~gsub("([A-z]+[0-9]+_[0-9]+)","\\1_bwa", .))
#
OUTPUT="./temp/data/asd.284_ped-gatk.txt"
fwrite(fam, file=OUTPUT, quote=FALSE, sep="\t", col.names = FALSE)
print(OUTPUT)