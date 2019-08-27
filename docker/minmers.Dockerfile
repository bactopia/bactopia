FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="0.0.5"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia-AP minmers"

COPY conda/minmers.yml /
RUN conda env create -f minmers.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-minmers/bin:$PATH
