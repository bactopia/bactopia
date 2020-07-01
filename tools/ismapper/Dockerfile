FROM nfcore/base

LABEL version="1.4.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image for Bactopia Tool ismapper"

COPY tools/ismapper/environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-ismapper/bin:$PATH
