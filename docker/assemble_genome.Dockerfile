FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="0.0.5"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia-AP assemble_genome"

COPY conda/assemble_genome.yml /
RUN conda env create -f assemble_genome.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-assemble_genome/bin:$PATH
