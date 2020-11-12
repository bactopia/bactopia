FROM nfcore/base

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia gather_fastqs process"

COPY conda/linux/gather_fastqs.yml /
RUN conda env create -f gather_fastqs.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-gather_fastqs/bin:$PATH
