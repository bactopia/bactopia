FROM nfcore/base

LABEL version="1.4.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia download_references process"

COPY conda/download_references.yml /
RUN conda env create -f download_references.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-download_references/bin:$PATH
