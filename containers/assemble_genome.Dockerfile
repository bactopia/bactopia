FROM nfcore/base:1.12.1

LABEL base.image="nfcore/base:1.12.1"
LABEL software="Bactopia - assemble_genome"
LABEL software.version="1.5.x"
LABEL description="A flexible pipeline for complete analysis of bacterial genomes"
LABEL website="https://bactopia.github.io/"
LABEL license="https://github.com/bactopia/bactopia/blob/master/LICENSE"
LABEL maintainer="Robert A. Petit III"
LABEL maintainer.email="robert.petit@emory.edu"

COPY conda/linux/assemble_genome.yml /
RUN conda env create -q -f assemble_genome.yml && conda clean -y -a
ENV PATH /opt/conda/envs/bactopia-assemble_genome/bin:$PATH
