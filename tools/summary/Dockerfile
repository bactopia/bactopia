FROM nfcore/base:1.12.1

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool summary"

COPY tools/summary/environment-linux.yml /
RUN conda env create -f /environment-linux.yml && conda clean -y -a
ENV PATH /opt/conda/envs/summary/bin:$PATH
