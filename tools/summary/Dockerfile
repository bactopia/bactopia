FROM nfcore/base

LABEL version="1.2.4"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool summary"

RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/summary/bin:$PATH
