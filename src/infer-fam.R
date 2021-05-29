#!/usr/bin/env Rscript
args=commandArgs(trailingOnly = TRUE)
#
library(readxl)
library(tidyverse)
# read info from excel
# INPUT1="./data/Sample-information.xlsx"
INPUT1=args[1]
df.info <- read_excel(INPUT1, range="A3:B104",
                      col_names = c("ID", "SEX"))
# insert info for ASD_014 due to missing sex info from excel file
df.info <- df.info %>%
  bind_rows(data.frame(ID="ASD_014_1", SEX="0"), .)
  
# extract FID
vt.fid <- gsub("_\\d$", "", df.info$ID) %>% unique() %>%
  gsub("ASD_", "ASD", x=.)
# mutate info into fam format knowing that the child has ASD and the parents are healthy
df.kid <- df.info %>%
  mutate(ID=gsub("ASD_", "ASD", x=ID)) %>%
  mutate(FID=gsub("_\\d$", "", ID)) %>%
  mutate(SEX=ifelse(SEX=="M",
                    1,
                    ifelse(SEX=="F", 2, 0))) %>%
  mutate(PHEN=2) %>%
  mutate(PID=gsub("_\\d$", "_3", ID)) %>%
  mutate(MID=gsub("_\\d$", "_2", ID)) %>%
  rename("IID"="ID") %>%
  select(c("FID", "IID", "PID", "MID","SEX", "PHEN"))
# add record of mom in FAM format
df.mom <- tibble(FID=vt.fid, IID=gsub("$", "_2", x=vt.fid), PID="0", MID="0", SEX=2, PHEN=1)
# add records of dad in FAM format
df.dad <- tibble(FID=vt.fid, IID=gsub("$", "_3", x=vt.fid), PID="0", MID="0", SEX=1, PHEN=1)
# combine records in FAM format
df.comb <- bind_rows(df.kid, df.mom, df.dad) %>%
  arrange(FID, IID)
# read  from fam
# INPUT2="./temp/data/asd.hg38_QC2.fam.bu"
INPUT2=args[2]
df.fam <- read_delim(INPUT2, delim="\t",
                     col_names = c("FID", "IID", "PID", "MID","SEX", "PHEN"))
# match id with df.comb
output <- df.fam %>%
  mutate(IID=paste(FID, IID, sep="_")) %>%
  select(1:2) %>%
  left_join(., df.comb[,-1], by="IID")
# create table in FAM format of plink
# OUTPUT="./temp/data/r-output.fam"
OUTPUT=args[3]
write_delim(output, file=OUTPUT, delim="\t",
            col_names = FALSE)
#
