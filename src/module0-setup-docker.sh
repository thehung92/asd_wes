#!/usr/bin/env Rscript
# get images
docker pull thehung92phuyen/biotools:v3.5
# run container at working directory
# cd [path/to/working directory]
# and mount the folder that have the reference of hg38.fasta and genetic map
# change REFDIR="/full/path/dir/contain-hg38.fasta-genetic.map/"
REFDIR="/Users/hung/Tools/Library/"
# REFDIR="/mnt/f/Data/Genomics"
docker run --name asd_wes \
    -v $PWD:/home/public/ \
    -v $REFDIR:/home/reference/ \
    -it thehung92phuyen/biotools:v3.5 bash

# docker run --name asd_wes -v ${PWD}:/home/public/ -it thehung92phuyen/biotools:v2.0 bash
# work inside container, the working directory is below
cd /home/public