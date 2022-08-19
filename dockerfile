FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

# install dependencies from pip3

RUN apt update && \
    apt install -y python3 && \
    apt install -y python-biopython \
        python3-pip \
        python3-pysam \
        wget \
        curl \
        bedtools \
        libxml2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        samtools \
        pigz \
        git \
        bc \
        r-base \
        unzip && \
    pip3 install biopython \
         pickle5 \
         pandas

ENV DEBIAN_FRONTEND=noninteractive

RUN R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('ggupset',dependencies=TRUE, repos='http://cran.rstudio.com/')"