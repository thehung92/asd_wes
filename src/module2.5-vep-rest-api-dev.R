#!/usr/bin/env Rscript
#
library(data.table)
library(tidyverse)
library(jsonlite)
library(httr)
# extract variants info from tdt analysis
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
# ensembls vep rest api
server <- "https://rest.ensembl.org"
ext <- "/vep/homo_sapiens/region/?CADD=1&canonical=1"
r.bonf <- POST(paste(server, ext, sep = ""),
          content_type("application/json"), accept("application/json"),
          body = body)
stop_for_status(r.bonf, task="VEP rest api query")
content(r.bonf) %>% toJSON() %>%
  fromJSON() %>%
  select(-"colocated_variants") -> df3
# unnest column with list of single element
df4 <- df3 %>%
  relocate(transcript_consequences, .after=last_col()) %>%
  unnest(!"transcript_consequences") %>%
  unnest("transcript_consequences", names_sep=".") %>%
  unnest(cols=starts_with("transcript_consequences") & !matches("consequence_terms"), keep_empty=TRUE) %>%
  rowwise() %>% mutate_at(vars(matches("consequence_terms")), toString) %>%
  filter(transcript_consequences.canonical=="1") %>%
  filter(!is.na(transcript_consequences.gene_symbol))
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
#### OUTPUT for supplementary ####
OUTPUT="output/tdt/asd.290_tdt-vep-rest-api.txt"
output <- tibble(head=c(head1, head2, head3))
write_delim(output, OUTPUT, delim="\t", col_names=FALSE)
write_delim(df4, OUTPUT, delim="\t", append=TRUE)

#### OUTPUT for plot-table ####
# required field, chrom, pos, ref, alt from df2
output <- df2 %>%
  left_join(., df1[,c("SNP", "A1", "OR", "P")], by=c("ID"="SNP")) %>% # annotate with A1, OR from df1
  mutate_at("POS", as.character) %>%
  rename("P-value"="P")
# transform df4 for merging
output <- df4 %>%
  separate(input, into=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"), sep=" ") %>%
  select("CHROM","POS","transcript_consequences.cadd_phred", "transcript_consequences.gene_symbol") %>%
  rename("CADD-score"=3, "SYMBOL"=4) %>%
  left_join(output, ., by=c("CHROM", "POS"))
# download latest SFARI database
URL <- "https://gene.sfari.org//wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=genes"
CMD <- paste("wget -P ./data/ --content-disposition", URL)
system(CMD)
# read in the newly download database
CMD <- "ls -t ./data/SFARI* | head -1" # this will get the latest file path
INPUT <- system(CMD, intern=TRUE)
db.sfari <- fread(file=INPUT)
output <- db.sfari[,c("gene-symbol", "gene-score")] %>%
  rename("SFARI-score"="gene-score") %>%
  left_join(output, ., by=c("SYMBOL"="gene-symbol"))
# format output
output %>%
  mutate_at(c("REF", "ALT"), ~ifelse(.==A1, paste0(.,"*"), .)) %>%
  select(-"A1") %>%
  unite("CHR:POS", 1:2, remove=TRUE, sep=":") %>%
  unite("REF/ALT",3:4, remove=TRUE, sep="/") %>%
  select(2,1,3:8) -> output
# official output is inside: OUTPUT="output/plot-table/table1.csv"
# write
OUTPUT="output/plot-table/table1_dev.csv"
heading <- paste("# SFARI database version:", gsub(".*(\\d+-\\d+-\\d+release).*", "\\1", INPUT)) %>%
  as.data.frame()
write_delim(heading, OUTPUT, delim="\t", col_names = FALSE)
write_delim(output, OUTPUT, delim="\t", append=TRUE)
