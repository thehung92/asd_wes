#!/usr/bin/env Rscript
#
library(data.table)
library(tidyverse)
#
URL <- "https://gene.sfari.org//wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=genes"
CMD <- paste("wget -P ./data/ --content-disposition", URL)
system(CMD)
# r <- download.file(URL, FILE, method="curl")
#
CMD <- "ls -t ./data/SFARI* | head -1" ; system(CMD, intern=TRUE)
