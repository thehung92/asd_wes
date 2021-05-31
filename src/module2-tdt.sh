#!/usr/bin/env bash
# populate the output structure
mkdir -p ./output/tdt
# use output from QC data module
INPUT="./output/asd.hg38_ibd-clean"
OUTPUT="./output/tdt/asd.290_ibd-clean"
plink --bfile $INPUT --tdt --adjust --out $OUTPUT
# convert to suitable format for plotting
OUTPUT2="./output/tdt/asd.ibd.tdt"
cat $OUTPUT.tdt | # read input
  awk 'BEGIN{OFS="\t";} {$1=$1; print}' | # change space delimited to tab delimited
  awk -v OFS="\t" '{gsub(/23/,"X", $1); print}' | # change 23 to x
  gsed '2,${s/^/chr/}' > asd.ibd.t  # add chr to the beginning > append to file