#!/usr/bin/env Rscript
# library
library(tidyverse)
# create list of parents from fam file
# this will be used to compute LD of unrelated individuals
df.fam <- read.delim("temp/data/asd.hg38_QC2.fam", header=FALSE)
colnames(df.fam) <- c("FID","IID","PID","MID","SEX","PHEN")
as_tibble(df.fam) %>%
  filter(grepl("_2|_3", IID)) -> output
OUTPUT="temp/ld-r2/asd.unrelated.fam"
system("mkdir -p temp/ld-r2")
write_delim(output, file=OUTPUT, delim="\t",
            col_names = FALSE)
# import tdt file to extract variants ID of adjacent variants from significant variants
df.tdt <- read.table("output/tdt/asd.290_ibd-clean.tdt", header=TRUE) %>%
  arrange(P)
# range of extraction
range <- 5e5
# loop through variants to compute ld-r2
for (i in 1:4) {
  # extract snp id
  sel.chr <- df.tdt[i, "CHR"]
  sel.pos <- df.tdt[i,"BP"]
  df.tdt[,1:3] %>%
    filter( CHR == sel.chr & BP >= sel.pos-range & BP <= sel.pos+range) %>%
    select(2) -> output
  OUTPUT=paste0("./temp/ld-r2/variant",i,"-region.txt")
  write_delim(output, file=OUTPUT, delim="\t", col_names = FALSE)
  # compute ld-r2 with plink
  INPUT="temp/data/asd.hg38_QC2"
  ARG1="temp/ld-r2/asd.unrelated.fam"
  ARG2=OUTPUT
  OUTPUT2=gsub(".txt","_ldr2", ARG2)
  CMD=paste("plink --bfile", INPUT, "--keep", ARG1, "--extract", ARG2,
            "--r2 square0 --out", OUTPUT2)
  system(CMD)
  # assign snpid to matrix
  mat.r2 <- read.delim(paste0(OUTPUT2, ".ld"), header=FALSE)
  vt.snp <- read.delim(OUTPUT, header=FALSE, col.names="SNP") %>% unlist()
  colnames(mat.r2) <- vt.snp
  rownames(mat.r2) <- vt.snp
  OUTPUT3=paste0("temp/ld-r2/matrix.r2.unrelated.variant",i,".txt")
  write_delim(output, file=OUTPUT3, delim="\t", col_names = FALSE)
}

# import matrix and list of variant
matrix.r2 <- read.delim("~/Data/Autism_vinmec_coop/run_2/LDscore/variant1.region.ld.ld", header=FALSE)
list.snp <- read.delim("/Users/hung/Data/Autism_vinmec_coop/run_2/Gviz_annotate_2/variant1.region.txt", col.names="SNP", header=FALSE)
# append colnames rownames and write file
colnames(matrix.r2) <- list.snp$SNP
rownames(matrix.r2) <- list.snp$SNP
write.table(matrix.r2, file="matrix.r2.variant1.txt",
            sep="\t", quote=FALSE)

for (i in 1:4) {
  # loop through 4 file
  # import matrix and list of variant
  input0=paste0("~/Data/Autism_vinmec_coop/run_2/LDscore/variant",i,".region.ldr2.unrelated.ld")
  input1=paste0("/Users/hung/Data/Autism_vinmec_coop/run_2/Gviz_annotate_2/variant",i,".region.txt")
  matrix.r2 <- read.delim(input0, header=FALSE)
  list.snp <- read.delim(input1, col.names="SNP", header=FALSE)
  # append colnames rownames and write file
  colnames(matrix.r2) <- list.snp$SNP
  rownames(matrix.r2) <- list.snp$SNP
  output=paste0("matrix.r2.unrelated.variant",i,".txt")
  write.table(matrix.r2, file=output,
              sep="\t", quote=FALSE)
}
