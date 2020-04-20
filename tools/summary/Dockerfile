FROM nfcore/base

LABEL version="1.3.1"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool summary"

COPY tools/summary/environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/summary/bin:$PATH
