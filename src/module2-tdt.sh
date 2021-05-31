#!/usr/bin/env bash
# populate the output structure
mkdir -p ./output/tdt ./output/plot-table
# use output from QC data module
INPUT="./output/asd.hg38_ibd-clean"
OUTPUT="./output/tdt/asd.290_ibd-clean"
plink --bfile $INPUT --tdt --adjust --out $OUTPUT
# run module 3 as Rscript

Rscript src/module3-tdt-analysis.R