#!/usr/bin/env Rscript
# lib
library(tidyverse)
library(gprofiler2)
library(openxlsx)
# import gene set
df0 <- read.csv(file="./output/plot-table/table1.csv")
df1 <- read.delim("output/plot-table/TADA_denovoonly.txt")
#
GeneList <- c(df1$gene.id, df0$SYMBOL)
#
dbList <- c("GO", "KEGG", "REAC", "WP","HP")
#
gostRslt2 <- gost(query = GeneList, 
                  organism = "hsapiens",
                  ordered_query = FALSE, 
                  multi_query = FALSE,
                  significant = TRUE,
                  exclude_iea = TRUE, 
                  measure_underrepresentation = FALSE,
                  evcodes = TRUE, 
                  user_threshold = 0.05, # change threshold to get more result
                  correction_method = "fdr", # false detection rate
                  domain_scope = "annotated",
                  custom_bg = NULL, 
                  numeric_ns = "",
                  sources = dbList,
                  as_short_link = FALSE)
#
as_tibble(gostRslt2$result) %>%
  select(c("source","term_id","term_name","term_size","intersection","p_value")) %>%
  group_by(source) %>%
  arrange(p_value, .by_group=TRUE) -> df
df %>% filter(p_value<=0.05) -> output
OUTPUT="output/plot-table/gprofiler_fdr.xlsx"
write.xlsx(output, file=OUTPUT)
