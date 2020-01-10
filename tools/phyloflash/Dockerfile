FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="1.2.4"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool phyloflash"

RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-phyloflash/bin:$PATH
