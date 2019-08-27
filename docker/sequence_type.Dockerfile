FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="0.0.5"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia-AP sequence_type"

COPY conda/sequence_type.yml /
RUN conda env create -f sequence_type.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-sequence_type/bin:$PATH
