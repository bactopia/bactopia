FROM nfcore/base

LABEL version="1.3.1"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia minmers process"

COPY conda/minmers.yml /
RUN conda env create -f minmers.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-minmers/bin:$PATH
