FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
LABEL authors="robert.petit@emory.edu" \
    description="Container image containing requirements for the Bactopia-AP insertion_sequences"

COPY conda/insertion_sequences.yml /
RUN conda env create -f insertion_sequences.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-insertion_sequences/bin:$PATH
