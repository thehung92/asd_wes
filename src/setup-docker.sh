#!/usr/bin/env Rscript
# save images from current container
docker commit -m 'first commit' db0c34a21bcd thehung92phuyen/biotools:v1.0 
# get images
docker pull thehung92phuyen/biotools:v1.0
# run container at working directory
# cd [path/to/working directory]
docker run --name asd_wes \
    -v $PWD:/home/public/ \
    -it thehung92phuyen/biotools:v1.0 bash
#
docker run --name asd_wes -v ${PWD}:/home/public/ -it thehung92phuyen/biotools:v1.0 bash
# work inside container