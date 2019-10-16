FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="1.2.0"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia gather_fastqs process"

COPY conda/gather_fastqs.yml /
RUN conda env create -f gather_fastqs.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-gather_fastqs/bin:$PATH
