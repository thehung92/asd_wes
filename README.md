# This project will analyze WES data from ASD cohort

## create data structure
```
# populate directory structure
mkdir -p data docker docs output/denovo output/plot-table output/tdt src temp/data temp/subset temp/ld-r2
```

## Instruction for running this project

* Set-up working environment by running the docker images or install the packages yourself
    * running docker images by running sr/module0-setup-docker.sh
    * the project is run inside the docker container
* running code from module1 to module7
    * the other scripts contain functions for each modules

## content of project dir after running on macOS
```
# tree .
├── data
├── docker
├── docs
├── output
│   ├── denovo
│   ├── plot-table
│   └── tdt
├── src
└── temp
    ├── data
    ├── ld-r2
    └── subset

12 directories
```
