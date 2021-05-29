#!/usr/bin/env Rscript
args=commandArgs(trailingOnly = TRUE)
#
library(data.table)
library(tidyverse)
#
#INPUT="F:/Data/Research/ASD_WES/temp/subset/asd.hg38_IBD.genome"
INPUT=args[1]
df0 <- fread(file=INPUT, header=TRUE)
# within family
df0 %>%
  select(1:4, "PI_HAT") %>%
  filter(FID1==FID2) -> df1
# within family and exlucde mom-dad pairs
df1 %>%
  filter(!(grepl("_2$", IID1) & grepl("_3$", IID2))) -> df2
# exclude accepted pi-hat number between 0.4-0.6
df2 %>%
  filter(PI_HAT<0.4 | PI_HAT>0.6) -> df3
# remove _4 _1 pair from list, it is identical twin
df3 %>%
  filter(!grepl("_4", IID2)) -> output
# OUTPUT=
OUTPUT=args[2]
fwrite(output, file=OUTPUT, sep="\t", quote=FALSE)