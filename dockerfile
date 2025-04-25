FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

# install dependencies from pip3

RUN apt update && \
    apt install -y python3 && \
    apt install -y \
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
         pandas \
         tqdm

ENV DEBIAN_FRONTEND=noninteractive

#RUN R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')"

#RUN R -e "install.packages('ggupset',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed && \
    chmod +x bigBedToBed && \
    mv bigBedToBed /usr/local/bin/bigBedToBed


# ARG UBUNTU_VER=18.04
# ARG CONDA_VER=latest
# ARG OS_TYPE=x86_64
# ARG PY_VER=3.9
# ARG PANDAS_VER=1.3

# # System packages 
# RUN apt-get update && apt-get install -yq curl wget jq vim

# # Use the above args 
# ARG CONDA_VER
# ARG OS_TYPE
# # Install miniconda to /miniconda
# RUN curl -LO "http://repo.continuum.io/miniconda/Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh"
# RUN bash Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh -p /miniconda -b
# RUN rm Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh
# ENV PATH=/miniconda/bin:${PATH}
# RUN conda update -y conda
# RUN conda init

# RUN conda install bioconda::ucsc-bigbedtobed
