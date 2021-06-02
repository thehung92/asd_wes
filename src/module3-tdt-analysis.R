#!/usr/bin/env Rscript
#
library(data.table)
library(tidyverse)
library(CMplot)
#
df0 <- fread(file="./output/tdt/asd.290_ibd-clean.tdt.adjusted")
df1 <- fread(file="./output/tdt/asd.290_ibd-clean.tdt")
# filter based on FDR_BH <0.05 with Region
as_tibble(df0) %>%
  filter(FDR_BH < 0.05) -> x
# annotate region with df1
x[,c(1,2,9)] %>%
  left_join(df1[,1:3], by=c("CHR", "SNP")) %>%
  select(1,4) %>%
  mutate(CHR=gsub("23", "X", CHR)) -> output
fwrite(x=output, file="output/tdt/regions.txt", sep="\t",
      quote=FALSE, col.names=FALSE)
print("output/tdt/regions.txt")
# extract variants info from vcf file
INPUT="./data/ASD_SAMPLES_bwa_gatk_variants_hg38_annotate.vcf.gz"
QUERY="bash ~/.bashrc ; bcftools query -R output/tdt/regions.txt -f '%CHROM %POS %ID %REF %ALT\n'"
CMD=paste(QUERY, INPUT)
df2 <- fread(cmd=CMD, header=FALSE, sep=" ",
             col.names = c("CHROM", "POS", "ID", "REF", "ALT"))
# write out in vcf format for annotation with VEP web tools
OUTPUT="./output/tdt/asd.290_tdt-sig-var.vcf"
df2 %>%
  mutate(QUAL=".", FILTER=".", INFO=".") %>%
  fwrite(x=., file=OUTPUT, quote=FALSE, sep="\t", col.names=FALSE)
cat(OUTPUT)
# annotate ./output/tdt/asd.290_tdt-sig-var.vcf with VEP web tools
# we will get asd.290_tdt-vep.txt
INPUT="./output/tdt/asd.290_tdt-vep.txt"

# take variants with BONF correction P-value 0.05 only
df0 %>%
  filter(BONF < 0.05) %>%
  select(1:3) %>%
  mutate(UNADJ=format(UNADJ, digits=3), CHR=gsub("23", "X", CHR)) -> output
# annotate output with chr:pos and ref-alt allele from vcf
output %>%
  left_join(., df2[,-1], by=c("SNP"="ID")) %>%
  select(2,1,4,5,6,3) -> output
# annotate output with A1 allele and OR from tdt
output %>%
  left_join(., df1[,c("SNP", "A1", "OR")], by="SNP") %>%
  mutate(OR=round(OR, digits=4)) -> output
# annotate output with gene and CADDscore from VEP
db.vep <- fread(file="./output/tdt/asd.290_tdt-vep.txt") %>% as_tibble()
db.vep[,c(1,6,37)] %>%
  rename(SNP=1) %>%
  left_join(output, ., by="SNP") %>%
  filter(!grepl("\\.", SYMBOL)) %>%
  mutate(CADD_PHRED=round(CADD_PHRED, 2)) -> output
# annotate with SFARI database
db.sfari <- fread(file="./data/SFARI-Gene_genes_10-29-2020release_11-04-2020export.csv")
db.sfari[,c("gene-symbol", "gene-score")] %>%
  left_join(output, ., by=c("SYMBOL"="gene-symbol")) -> output
# format output
output %>%
  mutate_at(c("REF", "ALT"), ~ifelse(.==A1, paste0(.,"*"), .)) %>%
  select(-"A1") %>%
  rename("P-value"=UNADJ, "SFARI-score"="gene-score", "CADD-score"="CADD_PHRED") %>%
  relocate("P-value", .after=last_col()) %>%
  unite("CHR:POS",2:3, remove=TRUE, sep=":") %>%
  unite("REF/ALT", 3:4, remove=TRUE, sep="/") -> output
# write output to table1.csv
OUTPUT="./output/plot-table/table1.csv"
write_csv(output, file=OUTPUT)

#### plot manhattan ####
rm(list=ls())
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
