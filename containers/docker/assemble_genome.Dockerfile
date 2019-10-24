FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>

LABEL version="1.2.2"
LABEL authors="robert.petit@emory.edu"
LABEL description="Container image containing requirements for the Bactopia assemble_genome process"

COPY conda/assemble_genome.yml /
RUN conda env create -f assemble_genome.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-assemble_genome/bin:$PATH
