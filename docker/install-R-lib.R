#!/usr/bin/env Rscript
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)
install.packages('data.table', dependencies=TRUE)
install.packages('tidyverse', dependencies=TRUE)