#!/usr/bin/env Rscript
#
library(data.table)
library(tidyverse)
library(CMplot)
#
# import data
IBD.tdt <- read.csv("output/tdt/asd.290_ibd-clean.tdt", sep="")
IBD.tdt.adj <- read.csv("output/tdt/asd.290_ibd-clean.tdt.adjusted", sep="")
# draw manhattan plot with CMplot (data format: SNP-CHR-BP-P1-P2)
as_tibble(IBD.tdt) %>%
  select(c(2,1,3,10)) %>%
  mutate(CHR=as.character(CHR)) %>%
  mutate_at("CHR", ~gsub(pattern="23", replacement="X", x=.)) -> df1
threshold <- 5e-2 / nrow(df1)
SNPs <- df1 %>% filter(P<=threshold) %>% select(1) %>% unlist()
# plot in R device
OUTPUT1="./output/plot-table/figure2.pdf"
pdf(file=OUTPUT1, width=8, height=4)
CMplot(df1,type="p",plot.type="m",LOG10=TRUE,chr.labels.angle=45,
       cex=0.5, amplify=TRUE, cex.lab=1, cex.axis=0.8,
       highlight=SNPs, highlight.text=SNPs, highlight.text.xadj=rep(0,4), highlight.text.cex=0.8,
       threshold=threshold,
       file="pdf",memo="", file.output=FALSE,verbose=TRUE,width=8,height=4)
dev.off()
# draw qq plot
OUTPUT2="./output/plot-table/figureS1.pdf"
pdf(file=OUTPUT2, width=4, height=4)
CMplot(df1,plot.type="q",box=FALSE, main=NULL,
       conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
       file="pdf",memo="",
       file.output=FALSE,verbose=TRUE,width=4,height=4)
dev.off()
