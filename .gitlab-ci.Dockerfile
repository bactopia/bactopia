FROM continuumio/miniconda3

LABEL version="1.4.10"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Gitlab CI"

RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*
RUN conda create -y -n bactopia -c conda-forge -c bioconda \
    ariba \
    beautifulsoup4 \
    biopython \
    blast \
    cd-hit \
    conda \
    coreutils \
    executor \
    lxml \
    mash \
    ncbi-genome-download \
    nextflow \
    pysam=0.15.3 \
    python=3.6.7 \
    requests \
    sed \
    unzip \
    wget \ 
    && conda clean -a && mkdir /opt/bactopia-envs

ENV PATH /opt/conda/envs/bactopia/bin:$PATH
