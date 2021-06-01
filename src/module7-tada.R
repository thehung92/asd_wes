#!/usr/bin/env Rscript
# lib
library(tidyverse)
library(readxl)
library(data.table)
source("src/TADA.v.1.2.R")
# function
vep_location <- function(chrom,pos,ref,alt) {
  #
  if (nchar(alt)-nchar(ref) < 0) {
    Location=paste0(chrom,":",pos,"-",as.numeric(pos)+nchar(ref)-nchar(alt))
  } else {
    Location=paste0(chrom,":",pos,"-",pos)
  }
  return(Location)
}
# 
# import mutation rate
df.dnmr <- read_excel("./data/gene_mutationrate.xls")
df.dnmr %>%
  mutate_at("splice_site", as.numeric) %>%
  mutate_at(c(4:9), function(x) {y <- 10^x}) %>%
  replace_na(list(splice_site=0)) -> df.dnmr
# compute mut.cls1=mis, mut.cls2=non+splice_site+frameshift
df.dnmr %>%
  select(c(2,6:9)) %>%
  mutate(mut.cls1=mis, mut.cls2=non+splice_site+frameshift) %>%
  select(c(1,6,7)) %>%
  rename(SYMBOL=1) -> df.dnmr

#### process denovo data ####
# import annotation
df.dn.vep <- read.delim("./output/denovo/asd.284_denovo-vep.txt")
df.dn.vep <- df.dn.vep %>%
  select("Location","Allele","Consequence","SYMBOL","CADD_PHRED") %>%
  mutate(Location=gsub(pattern="^",replacement="chr",x=Location)) %>%
  distinct() # remove duplicated row
# import vcf
INPUT="./output/denovo/asd.284_hiconfdenovo.vcf.gz"
CMD=paste("source ~/.bashrc; bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/hiConfDeNovo\t%INFO/AC\n'",INPUT)
df.dn.vcf <- fread(cmd=CMD, sep="\t", header=FALSE,
                   col.names=c("chrom","pos","ref","alt","denovo","AC")) 
# count number of mutation in vcf per sites by split with ","
df.dn.vcf$count <- apply(df.dn.vcf, 1, function(x) {
  x["denovo"] %>%
    strsplit(.,split=",") %>%
    unlist() %>% length()
}
)
# add Location to df
df.dn.vcf$Location <- apply(df.dn.vcf,1, function(x) {
  vep_location(x["chrom"], x["pos"], x["ref"], x["alt"])
}
) %>% gsub(pattern="\\s+",replacement="",x=.)
# merge the denovo count with vep output by "Location"
merge(df.dn.vcf, df.dn.vep, by="Location", all.x=TRUE) %>%
  as_tibble() -> df.dn.vcf.vep
#### classify as dmis
df.dn.vcf.vep %>%
  filter(grepl(pattern="missense", x=Consequence)) %>% # type missense
  filter(CADD_PHRED>=20) %>% # CADD >20
  filter(!grepl(pattern="splice", x=Consequence)) %>% # remove variants with splice consq
  as_tibble() -> dmis.dn
# count number of mutations per gene cls1=dmis
aggregate(dmis.dn$count, by=list(SYMBOL=dmis.dn$SYMBOL), FUN=sum) %>%
  rename(dn.cls1=2) -> dmis.dn
#### classify as LGD
df.dn.vcf.vep %>%
  filter(grepl(pattern="stop|frameshift|splice",x=Consequence)) %>%
  filter(!grepl(pattern="synonymous", x=Consequence)) %>%
  as_tibble() -> lgd.dn
# count number of mutations per gene cls2=lgd
aggregate(lgd.dn$count, by=list(SYMBOL=lgd.dn$SYMBOL), FUN=sum) %>%
  rename(dn.cls2=2) -> lgd.dn
#### finalize denovo input
# bind cols and replace NA with 0
merge(dmis.dn, lgd.dn, by="SYMBOL", all=TRUE) %>%
  replace_na(list(dn.cls1=0,dn.cls2=0)) -> tada.data
#### tada.data.1 ####

# merge with background mutation rate
merge(tada.data, df.dnmr,by="SYMBOL", all.x=TRUE) %>%
  as_tibble() -> tada.data
# remove NA in mut.class and rename gene.id; add case control column
# reformat tada.data
tada.data <- tada.data %>%
  rename(gene.id=1) %>%
  filter(!is.na(mut.cls1)) %>%
  mutate(case.cls1=0, ctrl.cls1=0, case.cls2=0, ctrl.cls2=0)
#### tada.data : ready ####

### specify the number of families included in the analysis
n.family = 96
n = data.frame(dn=n.family, ca=NA, cn=NA)
sample.counts <- list(cls1=n, cls2=n)

#### set up parameters####
lambda=2.0 # the burden
pi=0.05
gamma.mean.dn=(lambda-1)/pi+1
# calibrate
df.dnmr <- read_xls(path="~/Data/Autism_vinmec_coop/genemutation_framework/gene_mutationrate.xls")
# rate of denovo synonymous mutation across the genome (based on denovo mutation rate)
mu <- sum(10^df.dnmr$syn)
# expected number of synonysmous mutation
s.exp <- 2*n.family*mu
# observed number of denovo synonymous mutation, remove duplicated record due to vep annotation
df.dn.vep %>%
  filter(grepl(pattern="synonymous", x=Consequence)) %>%
  filter(!duplicated(Location)) %>% nrow() -> s.obs
# calibrate constant: s.obs/s.exp
s.obs/s.exp
# change mut.cls and mut.cls2 by multiply the original vector with the above constant
tada.data %>%
  mutate(mut.cls1=mut.cls1*(s.obs/s.exp)*0.32) %>%
  mutate(mut.cls2=mut.cls2*(s.obs/s.exp)) -> tada.data
#### edited in tada.data ####

# create the mutational data used by TADA-Denovo
cls1.counts=data.frame(dn=tada.data$dn.cls1, ca=NA, cn=NA)
rownames(cls1.counts)=tada.data$gene.id
cls2.counts=data.frame(dn=tada.data$dn.cls2, ca=NA, cn=NA)
rownames(cls2.counts)=tada.data$gene.id
tada.counts=list(cls1=cls1.counts,cls2=cls2.counts)

### set up mutation rates
mu=data.frame(cls1=tada.data$mut.cls1,cls2=tada.data$mut.cls2)

### set up denovo only TRUE/FALSE, here we WANT to restrict ourselves to de novo only analyses
denovo.only=data.frame(cls1=TRUE,cls2=TRUE)

# set up parameters
cls1= data.frame(gamma.mean.dn= 4.7,beta.dn=1,
                 gamma.mean.CC=NA,beta.CC=NA,rho1=NA,nu1=NA,rho0=NA,nu0=NA)
cls2= data.frame(gamma.mean.dn=20.0,beta.dn=1,
                 gamma.mean.CC=NA,beta.CC=NA ,rho1=NA,nu1=NA,rho0=NA,nu0=NA)
hyperpar=list(cls1=cls1,cls2=cls2)

# running TADA-Denovo
re.TADA <- do.call(cbind.data.frame, TADA(tada.counts=tada.counts, sample.counts=sample.counts, mu=mu, hyperpar=hyperpar, denovo.only=denovo.only))

# Bayesian FDR control
re.TADA$qval=Bayesian.FDR(re.TADA$BF.total, pi0 = 0.95)

# run permutation to get the null distributions to use for calculating p-values for TADA
re.TADA.null=do.call(cbind.data.frame, TADAnull(tada.counts=tada.counts, sample.counts=sample.counts, mu=mu, hyperpar=hyperpar, denovo.only=denovo.only, nrep=100))
re.TADA$pval=bayesFactor.pvalue(re.TADA$BF.total,re.TADA.null$BFnull.total)

# write out result table with qval < 0.05
re.TADA %>% arrange(qval) %>% filter(qval < 0.05) %>%
  mutate(gene.id=rownames(.)) %>%
  select(6,1:5) -> re.TADA.select
merge(re.TADA.select[,c(1,5,6)], tada.data[,1:3], by="gene.id", all.x=TRUE) %>%
  arrange(qval) %>%
  rename(DMIS=4, LGD=5) -> output
OUTPUT="./output/plot-table/table2.csv"
write_csv(output, file=OUTPUT)
#
#### read SFARI ####
asd.risk.gene <- read.csv("./data/SFARI-Gene_genes_10-29-2020release_11-04-2020export.csv")
as_tibble(asd.risk.gene) %>%
  select(2,6,7,9) %>%
  filter(gene.score > 0) -> asd.risk.gene
# write table S1
merge(tada.data[,1:3], asd.risk.gene, by.x="gene.id", by.y="gene.symbol") %>%
  rename(DMIS=2, LGD=3) -> output2
OUTPUT2="./output/plot-table/tables1.csv"
write_csv(output2, file=OUTPUT2)

