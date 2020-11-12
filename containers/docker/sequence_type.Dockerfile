FROM nfcore/base

LABEL version="1.5.x"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia sequence_type process"

COPY conda/linux/sequence_type.yml /
RUN conda env create -f sequence_type.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-sequence_type/bin:$PATH
