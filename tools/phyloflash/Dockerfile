FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="1.2.4"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool phyloflash"

RUN conda env create -f /environment.yml && conda clean -a \
    && apt-get -qq update \
    && apt-get -qq -y --no-install-recommends install python3 \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log /tmp/* /var/tmp/*
ENV PATH /opt/conda/envs/bactopia-phyloflash/bin:$PATH
