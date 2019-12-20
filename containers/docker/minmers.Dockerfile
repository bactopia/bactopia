FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="1.2.4"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia minmers process"

COPY conda/minmers.yml /
RUN conda env create -f minmers.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-minmers/bin:$PATH
