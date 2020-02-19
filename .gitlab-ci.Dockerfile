FROM continuumio/miniconda3

LABEL version="1.3.0"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Gitlab CI"

RUN conda create -y -n bactopia -c conda-forge -c bioconda \
    ariba \
    beautifulsoup4 \
    biopython \
    blast \
    cd-hit \
    conda \
    executor \
    lxml \
    mash \
    ncbi-genome-download \
    nextflow \
    pysam=0.15.3 \
    python=3.6.7 \
    requests \
    unzip \
    wget \ 
    && conda clean -a

ENV PATH /opt/conda/envs/bactopia/bin:$PATH
