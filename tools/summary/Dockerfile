FROM nfcore/base

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool summary"

COPY tools/summary/environment-linux.yml /
RUN conda env create -f /environment-linux.yml && conda clean -a
ENV PATH /opt/conda/envs/summary/bin:$PATH
