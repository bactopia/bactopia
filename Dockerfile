FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="1.0.1"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia"

RUN conda create -n bactopia -c conda-forge -c bioconda bactopia && conda clean -a
ENV PATH /opt/conda/envs/bactopia/bin:$PATH
