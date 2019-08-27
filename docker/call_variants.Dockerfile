FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="0.0.5"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia-AP call_variants"

COPY conda/call_variants.yml /
RUN conda env create -f call_variants.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-call_variants/bin:$PATH
