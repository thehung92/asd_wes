#!/usr/bin/env Rscript
#
library(data.table)
library(tidyverse)
library(jsonlite)
library(httr)
#
# extract variants info from tdt analysis
#
df0 <- fread(file="./output/tdt/asd.290_ibd-clean.tdt.adjusted")
df1 <- fread(file="./output/tdt/asd.290_ibd-clean.tdt")
# filter based on FDR_BH <0.05 with Region
as_tibble(df0) %>%
  filter(BONF < 0.05) -> df.bonf
# annotate region with df1
df.bonf[,c(1,2,9)] %>%
  left_join(df1[,1:3], by=c("CHR", "SNP")) %>%
  select(1,4) %>%
  mutate(CHR=gsub("23", "X", CHR)) %>%
  mutate(concat=paste(CHR,BP, sep=":")) %>%
  select(3) %>% unlist() %>% paste(., collapse = ",") -> ARG
# extract variants info from vcf file
INPUT="./data/ASD_SAMPLES_bwa_gatk_variants_hg38_annotate.vcf.gz"
QUERY="source ~/.bashrc ; bcftools query -f '%CHROM %POS %ID %REF %ALT\n' -r"
CMD=paste(QUERY, ARG, INPUT)
df2 <- fread(cmd=CMD, header=FALSE, sep=" ",
             col.names = c("CHROM", "POS", "ID", "REF", "ALT"))
# create vcf format for variants 
df2[,-3] %>%
  mutate(variants=paste(CHROM,POS,".",REF,ALT,".",".",".")) %>%
  select(5) %>% as.list() %>% toJSON() -> body
server <- "https://rest.ensembl.org"
ext <- "/vep/homo_sapiens/region/?CADD=1&canonical=1"
r.bonf <- POST(paste(server, ext, sep = ""),
          content_type("application/json"), accept("application/json"),
          body = body)
stop_for_status(r.bonf, task="VEP rest api query")
content(r.bonf) %>% toJSON() %>%
  fromJSON() %>%
  select(-"colocated_variants") -> df3
#
df4 <- apply(df3, 1, function(x) {
  # x is one row of df3 dataframe
  x["transcript_consequences"][[1]] %>%
    select(-"strand") %>%
    filter(canonical=="1") %>%
    filter(gene_symbol!="NULL") -> y2
  as_tibble(x) %>%
    select(-"transcript_consequences") %>%
    distinct() -> y1
  y <- cbind(y1, y2)
  return(y)
}) %>% bind_rows(.)
# save the version of rest api at the day of performance
server <- "https://rest.ensembl.org"
ext <- "/info/eg_version?"
v.info <- GET(paste(server, ext, sep = ""), content_type("application/json"))
toJSON(content(v.info)) %>%
  fromJSON()
