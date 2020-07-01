FROM nfcore/base

LABEL version="1.4.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia count_31mers process"

COPY conda/count_31mers.yml /
RUN conda env create -f count_31mers.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-count_31mers/bin:$PATH
