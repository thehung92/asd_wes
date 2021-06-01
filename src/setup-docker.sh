#!/usr/bin/env Rscript
# save images from current container
docker commit -m 'first commit' db0c34a21bcd thehung92phuyen/biotools:v1.0 
# get images
docker pull thehung92phuyen/biotools:v2.0
# run container at working directory
# cd [path/to/working directory]
# and mount the folder that have the reference of hg38 and genetic map
docker run --name asd_wes \
    -v $PWD:/home/public/ \
    -v /Users/hung/Tools/Library/:/home/reference/ \
    -it thehung92phuyen/biotools:v2.0 bash
#
docker run --name asd_wes -v ${PWD}:/home/public/ -it thehung92phuyen/biotools:v2.0 bash
# work inside container