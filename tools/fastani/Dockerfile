FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="1.2.2"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool ani"

RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/ani/bin:$PATH
