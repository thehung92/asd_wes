#!/usr/bin/env Rscript
#
sub_chr_ensembl <- function(x) {
  # x is a string or a vector of string
  x <- gsub("23","X", x)
  x <- gsub("24","Y", x)
  x <- gsub("25","MT",x)
  return(x)
}
#
sub_strand_gviz <- function(x) {
  #x is a vector of character
  x <- gsub("^1","+",x)
  x <- gsub("^-1","-",x)
  x <- gsub("0","*",x)
  #x <- factor(x, levels=c("+","-","*"))
  return(x)
}

