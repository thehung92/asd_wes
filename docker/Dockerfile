FROM continuumio/miniconda3:latest
LABEL MAINTAINER="thehung, thehung92phuyen@gmail.com"
# create virtual environment and install dependencies using conda
COPY . /root/tools/
WORKDIR /root/tools/
RUN conda create -y --name biotools
# make RUN commands use the new environment
SHELL ["conda", "run", "-n", "biotools", "/bin/bash", "-c"]
# RUN in the new environment to install R and library
RUN conda install -y -c conda-forge/label/main r-base \
    && Rscript install-R-lib.R
# install R library for genomic visualization, this require a lot of time, so a separate layer
RUN Rscript install-R-lib-2.R
# RUN in the new environment to install other biotools via conda
RUN conda install -y -c bioconda/label/main bcftools \
    && conda install -y -c bioconda/label/main plink \
    && conda install -y -c bioconda/label/main gatk4
# RUN to download and install plink2 directly from web
RUN apt-get update \
    && apt-get install -y wget \
    && apt-get install -y zip \
    && apt-get install -y less \
    && apt-get install -y nano \
    && wget https://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip \
    && unzip plink2_linux_x86_64.zip ; rm plink2_linux_x86_64.zip
# make bash start from biotools env and export /root/tools to $PATH, this layer define config
RUN conda init \
    && echo "conda activate biotools" >> ~/.bashrc \
    && echo "export PATH=/root/tools/:$PATH" >> ~/.bashrc