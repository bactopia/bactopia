FROM nfcore/base
MAINTAINER Robert A. Petit III <robert.petit@emory.edu>
LABEL authors="robert.petit@emory.edu" \
    description="Container image containing requirements for the Bactopia-AP annotate_genome"

COPY conda/annotate_genome.yml /
RUN conda env create -f annotate_genome.yml && conda clean -a
ENV PATH /opt/conda/envs/bactopia-annotate_genome/bin:$PATH
