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
* running code from module1 to module8 consecutively
    * the other scripts contain functions for each modules
    * .sh file should be run with bash shell
    * .R file should be run with R console
    * code should be run by each lines and following the comment inside each module

## content of project dir after running on macOS
```bash
# tree -d
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
