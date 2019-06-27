FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
LABEL authors="robert.petit@emory.edu" \
    description="Container image containing requirements for the Bactopia-AP download_references"

COPY conda/download_references.yml /
RUN conda env create -f download_references.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-download_references/bin:$PATH
