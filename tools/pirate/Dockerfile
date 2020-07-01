FROM nfcore/base

LABEL version="1.4.0"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool pirate"

COPY tools/pirate/environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-pirate/bin:$PATH
