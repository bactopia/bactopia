FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="0.0.5"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia-AP count_31mers"

COPY conda/count_31mers.yml /
RUN conda env create -f count_31mers.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-count_31mers/bin:$PATH
