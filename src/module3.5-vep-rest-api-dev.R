#!/usr/bin/env Rscript
#
library(data.table)
library(tidyverse)
library(jsonlite)
library(httr)
library(tidyjson)
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

#### test ####
# unnest column with list of single element
df3 %>%
  relocate(transcript_consequences, .after=last_col()) %>%
  unnest(!"transcript_consequences") %>%
  unnest("transcript_consequences", names_sep=".") %>% view()
  filter(transcript_consequences.canonical=="1") %>%
  filter(transcript_consequences.gene_symbol!="NULL") %>%
  unnest(cols=starts_with("transcript_consequences"), keep_empty=TRUE) -> temp
# str(temp, max.level=2)
# unnest nested column
df3 %>%
  select("transcript_consequences") %>%
  unnest(everything(), keep_empty=TRUE)
# tidyjson
content(r.bonf) -> temp
as_tbl_json(temp) %>%
#### ####
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
# convert df4 list value to write to file
df4 %>%
  mutate_if(is.list, function(x){
    unlist(x) %>% paste(., collapse=",")
  }) -> temp
# save the version of rest api at the day of running code and write to output
server <- "https://rest.ensembl.org"
ext <- "/info/eg_version?"
v.info <- GET(paste(server, ext, sep = ""), content_type("application/json"))
stop_for_status(v.info, task="VEP rest api query")
toJSON(content(v.info)) %>%
  fromJSON() -> temp
head1 <- paste("# Ensembl Genomes version of the database:", names(temp), temp[[1]])
#
ext <- "/info/rest?"
v.info <- GET(paste(server, ext, sep = ""), content_type("application/json"))
stop_for_status(v.info, task="VEP rest api query")
toJSON(content(v.info)) %>%
  fromJSON() -> temp
head2 <- paste("# Version of the Ensembl REST API:", names(temp), temp[[1]])
#
ext <- "/info/software?"
v.info <- GET(paste(server, ext, sep = ""), content_type("application/json"))
stop_for_status(v.info, task="VEP rest api query")
toJSON(content(v.info)) %>%
  fromJSON() -> temp
head3 <- paste("# Version of the Ensembl API used by the REST server:", names(temp), temp[[1]])
# OUTPUT
OUTPUT="output/tdt/asd.290_tdt-vep-rest-api.txt"
output <- tibble(head=c(head1, head2, head3))
write_delim(output, OUTPUT, delim="\t", col_names=FALSE)
write_delim(df4, OUTPUT, delim="\t", append=TRUE)