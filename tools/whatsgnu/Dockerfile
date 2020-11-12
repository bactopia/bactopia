FROM nfcore/base

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool whatsgnu"

COPY tools/whatsgnu/environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-whatsgnu/bin:$PATH
