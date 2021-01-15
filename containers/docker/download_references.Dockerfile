FROM nfcore/base:1.12.1

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia download_references process"

COPY conda/linux/download_references.yml /
RUN conda env create -f -q download_references.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-download_references/bin:$PATH
